import asyncio
from typing import (
    Awaitable,
    Callable,
    List,
    Dict
)

import aiohttp

from build_module.chembl_data_processing.services.api_client import ChemblAPIClient
from build_module.chembl_data_processing.services.db_handler import DatabaseHandler
from build_module.fingerprint_computation.fingerprint_processor import FingerprintProcessor
from build_module.utils.config import (
    MOLECULE_URL,
    CHEMBL_ID_LOOKUP_URL,
)
from build_module.utils.logging_config import setup_logging

logger = setup_logging()


async def process_and_insert_data(
        api_client: ChemblAPIClient,
        process_func: Callable[[List[Dict]], Awaitable[None]]
) -> None:
    """
    Fetch data using the provided API client and process it using the given function.

    Args:
        api_client (ChemblAPIClient): The API client to use for fetching data.
        process_func (Callable): The function to process and insert the fetched data.
    """
    async with aiohttp.ClientSession() as session:
        try:
            data_chunks = await api_client.fetch_data_chunks(session)
            for chunk in data_chunks:
                await process_func(chunk)
        except Exception as e:
            logger.error(f'Error in process_and_insert_data: {e}')


async def main() -> None:
    """
    Main function to orchestrate the ETL process.
    """
    db_handler = DatabaseHandler()
    await db_handler.setup()

    chembl_id_client = ChemblAPIClient(CHEMBL_ID_LOOKUP_URL)
    try:
        await process_and_insert_data(
            chembl_id_client,
            db_handler.process_chembl_id_lookup_data,
        )
    except Exception as em:
        logger.error(f'Error happened during processing staging_chembl_id_lookup data : {em}')

    molecule_client = ChemblAPIClient(MOLECULE_URL)
    try:
        await process_and_insert_data(
            molecule_client,
            db_handler.process_molecule_data,
        )
    except Exception as en:
        logger.error(f'Error happened during processing molecules : {en}')


def fingerprint_computation():
    """
    Here we compute fingerprints based on database staging_compound_structures data
    """
    fingerprint_processor = FingerprintProcessor()
    fingerprint_processor.process_and_store_fingerprints(table_name='staging_compound_structures')


if __name__ == '__main__':
    try:
        asyncio.run(main())
        fingerprint_computation()
    except Exception as e:
        logger.error(f'An error occurred during the ETL process: {e}')
