import re
from typing import List

import pandas as pd

from utils.logging_config import setup_logging
from utils.s3_operations import AsyncS3Saver

logger = setup_logging()


class ProcessTargetMolecules:
    def __init__(
            self,
            s3_saver: AsyncS3Saver,
            input_folder: str,
    ):
        self.s3_saver = s3_saver
        self.input_folder = input_folder

    def load_and_clean_csv_files(
            self,
            file_list: List[str],
    ) -> pd.DataFrame:

        logger.info(f'Processing {len(file_list)} files from S3 bucket')

        dfs: List[pd.DataFrame] = []
        for csv_file in file_list:
            try:
                decoded_content = self.s3_saver.read_csv_from_s3(csv_file)
                if decoded_content.strip():
                    cleaned_df = self._clean_csv(
                        decoded_content,
                        csv_file,
                    )
                    dfs.append(cleaned_df)
                else:
                    logger.warning(f'Empty file detected: {csv_file}')
            except Exception as e:
                logger.error(f'Error processing file {csv_file}: {str(e)}')

        mols_df = pd.concat(
            dfs,
            ignore_index=True,
        )

        duplicates = mols_df[mols_df.duplicated(
            subset='chembl_id',
            keep=False,
        )]
        for _, row in duplicates.iterrows():
            logger.warning(f'Duplicate ChEMBL ID found. Keeping only the first occurrence.')

            mols_df.drop_duplicates(
                subset='chembl_id',
                keep='first',
                inplace=True,
            )
        logger.info(f"Total number of molecules loaded and cleaned: {len(mols_df)}")

        return mols_df

    @staticmethod
    def _clean_csv(content: str, filename: str) -> pd.DataFrame:
        lines = content.split('\n')
        processed_lines = []

        for i, line in enumerate(lines[1:], start=2):
            if not line.strip():
                continue

            parts = line.split(',', 1)
            if len(parts) != 2:
                logger.warning(f'File {filename}, Line {i}: Skipping malformed line: {line}')
                continue

            chembl_id, smiles = parts[0].strip(), parts[1].strip()

            if re.search(r'[^\x00-\x7F]', line):
                logger.warning(f'File {filename}, Line {i}: Skipping line with unknown characters: {line}')
                continue

            if smiles.count(',') >= 1:
                logger.warning(f'File {filename}, Line {i}: Skipping line with multiple commas in SMILES: {line}')
                continue

            processed_lines.append((chembl_id, smiles))

        return pd.DataFrame(processed_lines, columns=['chembl_id', 'smiles'])
