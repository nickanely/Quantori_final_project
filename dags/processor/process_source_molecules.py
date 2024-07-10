import pandas as pd
import pyarrow.dataset as ds
from rdkit import DataStructs
from tenacity import retry, stop_after_attempt, wait_fixed

from utils.logging_config import setup_logging
from utils.s3_operations import AsyncS3Saver

logger = setup_logging()


class ProcessSourceMolecules:

    def __init__(
            self,
            s3_saver: AsyncS3Saver,
            s3_fingerprints_folder: str
    ):
        self.s3_saver = s3_saver
        self.s3_fingerprints_folder = s3_fingerprints_folder
        self.result_df = pd.DataFrame(columns=['chembl_id', 'computed_fingerprint'])

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_fixed(10)
    )
    def load_and_precompute_fingerprints(self) -> pd.DataFrame:

        fingerprint_files = self.s3_saver.list_files(
            self.s3_fingerprints_folder,
            file_extension='.parquet'
        )

        for file in fingerprint_files:
            logger.info(f'Processing file: {file}', )
            try:
                self._process_fingerprint_file(file)
            except Exception as exc:
                logger.error(f'Error processing file {file}: {exc}')

        self.result_df.drop_duplicates(
            subset='chembl_id',
            inplace=True
        )
        logger.info(f'Processed a total of {len(self.result_df)} unique fingerprints')
        return self.result_df

    @retry(stop=stop_after_attempt(5), wait=wait_fixed(10))
    def _process_fingerprint_file(
            self,
            file_path: str
    ) -> None:

        try:
            dataset = ds.dataset(
                file_path,
                format='parquet',
                filesystem=self.s3_saver.s3fs_client
            )
            for batch in dataset.to_batches(columns=['chembl_id', 'fingerprint']):
                df = batch.to_pandas()
                df['computed_fingerprint'] = df['fingerprint'].apply(self._create_bit_vect)
                chunk_result = df[
                    ['chembl_id', 'computed_fingerprint']
                ].dropna()
                self._append_chunk(chunk_result)

            logger.info(f'Finished processing fingerprint file: {file_path}')
        except Exception as e:
            logger.error(f'Failed to process fingerprint file {file_path}: {str(e)}')

    def _append_chunk(self, chunk_df: pd.DataFrame) -> None:
        self.result_df = pd.concat([self.result_df, chunk_df], ignore_index=True)

    @staticmethod
    def _create_bit_vect(fingerprint):
        if fingerprint is None:
            logger.warning('Received None as fingerprint, returning None.')
            return None
        try:
            bit_vect = DataStructs.cDataStructs.CreateFromBitString(fingerprint)
            logger.debug('Successfully created bit vector.')
            return DataStructs.BitVectToText(bit_vect)
        except Exception as e:
            logger.warning(f'Failed to create bit vector from fingerprint: {fingerprint[:50]}... Error: {str(e)}')
            return None
