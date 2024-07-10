import asyncio
import io
from typing import List, Tuple

import aioboto3
import pandas as pd
import s3fs

from  config import (
    S3_FINGERPRINTS_FOLDER_PATH,
    S3_SIMILARITY_FOLDER_PATH,
    S3_TOP10_FOLDER_PATH,
    S3_INPUT_FILES,
    S3_BUCKET_NAME,
    AWS_ACCESS_KEY,
    AWS_SECRET_KEY,
    AWS_SESSION_TOKEN,

)
from logging_config import setup_logging

logger = setup_logging()


class AsyncS3Saver:
    def __init__(self, max_workers: int = 10):

        self.s3_bucket_name = S3_BUCKET_NAME
        self.s3_top10_folder_path = S3_TOP10_FOLDER_PATH
        self.s3_similarity_folder_path = S3_SIMILARITY_FOLDER_PATH
        self.s3_fingerprints_folder_path = S3_FINGERPRINTS_FOLDER_PATH
        self.input_folder = S3_INPUT_FILES
        self.max_workers = max_workers
        self.s3fs_client = s3fs.S3FileSystem(
            key=AWS_ACCESS_KEY,
            secret=AWS_SECRET_KEY,
            token=AWS_SESSION_TOKEN,
        )
        self.session = aioboto3.Session(
            aws_access_key_id=AWS_ACCESS_KEY,
            aws_secret_access_key=AWS_SECRET_KEY,
            aws_session_token=AWS_SESSION_TOKEN,
        )
        self.chunk_counter = 0

    async def save_to_s3(self, buffer: io.BytesIO, key: str):
        try:
            async with self.session.client('s3') as s3:
                await s3.upload_fileobj(
                    buffer,
                    self.s3_bucket_name,
                    key,
                )
            logger.info(f'Successfully uploaded {key} to bucket {self.s3_bucket_name}.')
        except Exception as e:
            logger.error(f'Failed to upload {key} to bucket {self.s3_bucket_name}: {str(e)}')

    async def save_full_similarity_result(self, chembl_id: str, full_similarity_scores: pd.DataFrame):
        try:
            full_similarity_buffer = io.BytesIO()
            full_similarity_scores[
                ['chembl_id', 'similarity_score']
            ].to_parquet(
                full_similarity_buffer,
                index=False,
            )
            full_similarity_buffer.seek(0)
            full_similarity_key = f'{self.s3_similarity_folder_path}{chembl_id}_similarity_score.parquet'
            await self.save_to_s3(full_similarity_buffer, full_similarity_key)
        except Exception as e:
            logger.error(f'Failed to process full similarity result for {chembl_id}: {str(e)}')

    async def save_similarity_results(self, similarity_results: List[Tuple[str, pd.DataFrame, pd.DataFrame]]):
        tasks = []
        for chembl_id, top_10_similar, full_similarity_scores in similarity_results:
            tasks.append(
                self.save_full_similarity_result(
                    chembl_id,
                    full_similarity_scores,
                ))

        await asyncio.gather(*tasks)
        logger.info(f'Saved similarity results for {len(similarity_results)} molecules.')

    def save_similarity_results_in_parallel(self, similarity_results: List[Tuple[str, pd.DataFrame, pd.DataFrame]]):
        asyncio.run(self.save_similarity_results(similarity_results))

    def list_files(self, folder: str, file_extension: str = '') -> List[str]:
        try:
            files = self.s3fs_client.glob(
                f"{self.s3_bucket_name}/{folder}*{file_extension}"
            )
            logger.info(f'Found {len(files)} files with extension {file_extension} in {folder}')
            return files
        except Exception as e:
            logger.error(f'Error listing files in {folder}: {str(e)}')
            return []

    def read_csv_from_s3(self, s3_path: str) -> str:
        try:
            with self.s3fs_client.open(s3_path, 'rb') as f:
                content = f.read()
                try:
                    decoded_content = content.decode('utf-8')
                except UnicodeDecodeError:
                    try:
                        decoded_content = content.decode('latin-1')
                        logger.warning(f'File {s3_path} is not in UTF-8 encoding. Processed with latin-1 encoding.')
                    except UnicodeDecodeError:
                        logger.error(f'File {s3_path} has an unsupported encoding. Skipping this file.')

                return decoded_content
        except Exception as e:
            logger.error(f'Error reading CSV file from {s3_path}: {str(e)}')
            return ''
