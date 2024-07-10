import asyncio
from typing import (
    Optional,
    List,
    Dict,
)

import asyncpg
import pandas as pd

from build_module.chembl_data_processing.process_data import (
    process_chembl_id_lookup_chunk,
    process_molecule_chunk,
)
from build_module.utils.config import (
    DATABASE_URL,
    API_CHUNK_SIZE,
)
from build_module.utils.logging_config import setup_logging

logger = setup_logging()


class DatabaseHandler:
    """
    DatabaseHandler manages database interactions, including connection pooling and asynchronous data insertion.

    Attributes:
        db_url (str): Database connection URL.
        chunk_size (int): Size of data chunks for batch insertion.
        pool (asyncpg.pool.Pool): Connection pool for the database.
    """

    def __init__(self):
        self.db_url = DATABASE_URL
        self.chunk_size = API_CHUNK_SIZE
        self.pool: Optional[asyncpg.pool.Pool] = None

    async def setup(self) -> None:
        """
        Set up the database connection pool.
        """
        self.pool = await asyncpg.create_pool(self.db_url)

    async def insert_dataframe_to_db(self, df: pd.DataFrame, table_name: str) -> None:
        """
        Insert a DataFrame into the database in chunks asynchronously.

        Args:
            df (pd.DataFrame): DataFrame to insert into the database.
            table_name (str): Name of the database table to insert data into.
        """
        chunks = [df.iloc[i:i + self.chunk_size] for i in range(0, len(df), self.chunk_size)]
        try:
            await asyncio.gather(*[self.insert_chunk(chunk, table_name) for chunk in chunks])
        except Exception as e:
            logger.error(f'Error in async data insertion: {e}')

    async def insert_chunk(self, chunk_df: pd.DataFrame, table_name: str) -> None:
        """
        Insert a chunk of data into the database.

        Args:
            chunk_df (pd.DataFrame): DataFrame chunk to insert into the database.
            table_name (str): Name of the database table to insert data into.
        """
        try:
            chunk_df = chunk_df.where(pd.notnull(chunk_df), None)
            records = chunk_df.to_dict(orient='records')
            values = [tuple(record.values()) for record in records]
            columns = ', '.join(records[0].keys())
            placeholders = ', '.join([f'${i + 1}' for i in range(len(records[0]))])
            query = f'INSERT INTO {table_name} ({columns}) VALUES ({placeholders})'

            async with self.pool.acquire() as connection:
                async with connection.transaction():
                    await connection.executemany(query, values)

            logger.info(f'Inserted chunk into {table_name} - {len(chunk_df)} records')
        except Exception as e:
            logger.error(f'Error inserting data into {table_name}: {e}')

    async def process_chembl_id_lookup_data(self, aggregated_records: List[Dict]) -> None:
        """
        Process and insert ChEMBL ID lookup data into the database.
        """
        df = process_chembl_id_lookup_chunk(aggregated_records)
        await self.insert_dataframe_to_db(df, 'staging_chembl_id_lookup')

    async def process_molecule_data(self, aggregated_records):
        (molecule_dict_df,
         compound_properties_df,
         compound_structures_df) = process_molecule_chunk(aggregated_records)

        await asyncio.gather(
            self.insert_dataframe_to_db(molecule_dict_df, 'staging_molecule_dictionary'),
            self.insert_dataframe_to_db(compound_properties_df, 'staging_compound_properties'),
            self.insert_dataframe_to_db(compound_structures_df, 'staging_compound_structures')
        )
