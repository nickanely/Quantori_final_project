import concurrent.futures
from typing import Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from tqdm import tqdm

from config import (
    COMPUTE_NUM_PROCESSES,
    FPS_MOL_RADIUS,
    FPS_BITS
)
from utils.logging_config import setup_logging
from utils.s3_operations import AsyncS3Saver

logger = setup_logging()


class ProcessSimilarityScore:
    def __init__(
            self,
            s3_saver: AsyncS3Saver,
            target_mols_df: pd.DataFrame,
            source_mols_df: pd.DataFrame
    ):
        self.s3_saver = s3_saver
        self.target_mols_df = target_mols_df
        self.source_mols_df = source_mols_df
        self.source_fps = self.source_mols_df['computed_fingerprint'].tolist()

    @staticmethod
    def _parse_smiles(smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f'Invalid SMILES string: {smiles}')
        return mol

    def _calculate_similarity_score(self, target_fps) -> pd.DataFrame:
        similarity_scores = DataStructs.BulkTanimotoSimilarity(
            target_fps,
            self.source_fps
        )
        result_df = self.source_mols_df.copy()
        result_df['similarity_score'] = similarity_scores
        return result_df[
            ['chembl_id', 'similarity_score']
        ].sort_values('similarity_score', ascending=False)

    def _process_single_target(
            self,
            target_row: pd.Series
    ) -> Tuple[str, pd.DataFrame, pd.DataFrame]:

        target_smiles = target_row.smiles
        chembl_id = target_row.chembl_id
        logger.info(f'Processing target molecule: {chembl_id}')

        try:
            target_mol = self._parse_smiles(target_smiles)
            target_fps = AllChem.GetMorganFingerprintAsBitVect(
                target_mol,
                FPS_MOL_RADIUS, FPS_BITS
            )

            similarity_scores_df = self._calculate_similarity_score(target_fps)

            if similarity_scores_df.empty:
                logger.warning(f'No similarity scores calculated for {chembl_id}')
                return (
                    chembl_id,
                    pd.DataFrame(),
                    pd.DataFrame()
                )

            similarity_scores_df = similarity_scores_df[similarity_scores_df['chembl_id'] != chembl_id]

            top_10_similar = self._compute_top_10(similarity_scores_df)

            logger.info(f'Finished processing target molecule: {chembl_id}')

            return (
                chembl_id,
                top_10_similar,
                similarity_scores_df
            )

        except Exception as e:
            logger.error(f'Error processing target molecule {chembl_id}: {str(e)}')
            return (
                chembl_id,
                pd.DataFrame(),
                pd.DataFrame()
            )

    @staticmethod
    def _compute_top_10(similarity_scores_df):
        top_10_similar = similarity_scores_df.head(10).copy()
        top_10_similar['has_duplicates_of_last_largest_score'] = False

        if len(top_10_similar) == 10:
            last_score = top_10_similar.iloc[-1]['similarity_score']
            if (similarity_scores_df.iloc[10:]['similarity_score'] == last_score).any():
                top_10_similar.loc[
                    top_10_similar.index[-1],
                    'has_duplicates_of_last_largest_score'
                ] = True

        return top_10_similar

    def compute_similarity_scores(self):
        logger.info('Starting to compute similarity scores')

        chunksize = max(
            1,
            len(self.target_mols_df) // (COMPUTE_NUM_PROCESSES * 4)
        )

        with concurrent.futures.ThreadPoolExecutor(max_workers=COMPUTE_NUM_PROCESSES) as executor:
            futures = executor.map(
                self._process_single_target,
                self.target_mols_df.itertuples(index=False),
                chunksize=chunksize,
            )

            similarity_results = [
                result for result in tqdm(
                    futures,
                    total=len(self.target_mols_df),
                    desc='Processing',
                    unit='molecule'
                )
                if result[1] is not None and not result[1].empty and not result[2].empty
            ]

        logger.info(f'Similarity scores computed for {len(similarity_results)} target molecules')
        logger.info(f'Started uploading similarity scores to S3...')

        self.s3_saver.save_similarity_results_in_parallel(similarity_results)

        return [(result[0], result[1]) for result in similarity_results]
