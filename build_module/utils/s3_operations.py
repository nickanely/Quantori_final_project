import asyncio
import io
import logging
from typing import List, Tuple

import aioboto3
import botocore
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from build_module.utils.config import (
    AWS_SESSION_TOKEN,
    AWS_SECRET_KEY,
    AWS_ACCESS_KEY,
    S3_BUCKET_NAME,
    S3_FINGERPRINTS_FOLDER_PATH,
    COMPUTE_NUM_PROCESSES,
)


class AsyncS3Saver:
    def __init__(self):
        self.s3_bucket_name = S3_BUCKET_NAME
        self.s3_fingerprints_folder_path = S3_FINGERPRINTS_FOLDER_PATH
        self.max_workers = COMPUTE_NUM_PROCESSES
        self.session = aioboto3.Session(
            aws_access_key_id=AWS_ACCESS_KEY,
            aws_secret_access_key=AWS_SECRET_KEY,
            aws_session_token=AWS_SESSION_TOKEN,
        )

    async def save_to_s3(self, buffer: io.BytesIO, key: str):
        try:
            logging.info(f"Starting upload of {key} to bucket {self.s3_bucket_name}")
            config = botocore.config.Config(connect_timeout=300, read_timeout=300)
            async with self.session.client('s3', config=config) as s3:
                await s3.upload_fileobj(buffer, self.s3_bucket_name, key)
            logging.info(f"Successfully uploaded {key} to bucket {self.s3_bucket_name}.")
        except Exception as e:
            logging.error(f"Failed to upload {key} to bucket {self.s3_bucket_name}: {str(e)}")

    async def save_fingerprints(self, fingerprints: List[Tuple[int, pd.DataFrame]]):
        tasks = []
        counter = 0
        for chunk_id, df in fingerprints:
            buffer = io.BytesIO()
            table = pa.Table.from_pandas(df)
            pq.write_table(table, buffer)
            buffer.seek(0)
            key = f"{self.s3_fingerprints_folder_path}fingerprints_chunk_{counter}.parquet"
            tasks.append(self.save_to_s3(buffer, key))
            counter += 1

        await asyncio.gather(*tasks)
        logging.info(f"Saved fingerprints for {len(fingerprints)} chunks.")

    def save_fingerprints_in_parallel(self, fingerprints: List[Tuple[int, pd.DataFrame]]):
        asyncio.run(self.save_fingerprints(fingerprints))
        logging.info(f"Initiated saving fingerprints for {len(fingerprints)} chunks.")
