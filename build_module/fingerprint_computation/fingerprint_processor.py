import concurrent.futures
import math

from typing import Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import pandas as pd
from sqlalchemy import create_engine, text

from build_module.utils.s3_operations import AsyncS3Saver

from build_module.utils.config import (
    DATABASE_URL, FINGERPRINT_CHUNK_SIZE, COMPUTE_NUM_PROCESSES, FPS_MOL_RADIUS, FPS_BITS
)

from build_module.utils.logging_config import setup_logging

logger = setup_logging()


def compute_fingerprint(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f'Invalid SMILES: {smiles}')
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            FPS_MOL_RADIUS,
            FPS_BITS
        )
        return fp.ToBitString()
    except Exception as e:
        logger.error(f'Error computing fingerprint for SMILES {smiles}: {str(e)}')
        return None


def process_chunk(df: pd.DataFrame) -> pd.DataFrame:
    df['fingerprint'] = df['canonical_smiles'].apply(compute_fingerprint)
    return df[['chembl_id', 'fingerprint']].dropna(subset=['fingerprint'])


def fetch_chunk(engine, table_name: str, offset: int, limit: int) -> pd.DataFrame:
    query = text(f"""
        SELECT chembl_id, canonical_smiles
        FROM {table_name}
        LIMIT :limit OFFSET :offset
    """)
    with engine.connect() as connection:
        result = connection.execute(
            query,
            {'limit': limit, 'offset': offset},
        )
        return pd.DataFrame(
            result.fetchall(),
            columns=result.keys(),
        )


def fetch_and_process_chunk(args: Tuple[str, str, int, int]) -> Tuple[int, pd.DataFrame]:
    connection_string, table_name, offset, limit = args
    engine = create_engine(connection_string)

    with engine.connect() as connection:
        df = fetch_chunk(
            connection,
            table_name,
            offset,
            limit,
        )
        logger.info(f'Fetched chunk with offset {offset} and {len(df)} rows')

        if not df.empty:
            processed_df = process_chunk(df)

            return offset, processed_df

    return offset, pd.DataFrame()


class FingerprintProcessor:

    def __init__(self):
        self.engine = create_engine(DATABASE_URL)
        self.s3_saver = AsyncS3Saver()

    def get_total_rows(self, table_name: str) -> int:
        with self.engine.connect() as connection:
            count_query = text(f"SELECT COUNT(DISTINCT chembl_id) FROM {table_name}")
            return connection.execute(count_query).scalar()

    def process_and_store_fingerprints(self, table_name: str):
        total_rows = self.get_total_rows(table_name)
        logger.info(f'Total unique molecules in database: {total_rows}')

        num_chunks = math.ceil(total_rows / FINGERPRINT_CHUNK_SIZE)
        chunk_size = math.ceil(total_rows / num_chunks)

        chunk_args = [
            (DATABASE_URL,
             table_name,
             i * chunk_size,
             chunk_size,
             )
            for i in range(num_chunks)
        ]
        progress_bar = tqdm(
            total=num_chunks,
            desc='Processing',
            unit='chunk',
        )
        results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=COMPUTE_NUM_PROCESSES) as executor:
            futures = [executor.submit(
                fetch_and_process_chunk,
                arg,
            ) for arg in chunk_args]

            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results.append(result)
                logger.info(f'Processed chunk {len(results)} of {num_chunks}')
                progress_bar.update(1)

        progress_bar.close()

        total_processed = sum(len(df) for _, df in results)
        logger.info(f'Completed processing. Total fingerprints computed: {total_processed}')
        logger.info('started saving results in s3 bucket...')
        self.s3_saver.save_fingerprints_in_parallel(results)

        if total_processed != total_rows:
            logger.warning(
                f'Mismatch between unique molecules ({total_rows}) and processed fingerprints ({total_processed})')
