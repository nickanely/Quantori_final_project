import asyncio
import json
import time
from typing import (
    List,
    Dict,
    Optional
)

import aiohttp
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential
)

from build_module.utils.config import (
    BATCH_SIZE,
    TIMEOUT,
    MAX_RETRIES,
    MAX_CONCURRENT_REQUESTS,
    API_RATE_LIMIT,
)
from build_module.utils.logging_config import setup_logging

logger = setup_logging()


class RateLimiter:
    """
      RateLimiter ensures that API calls do not exceed the allowed rate.

      Attributes:
          calls_per_second (float): Maximum number of API calls per second.
          last_call (float): Timestamp of the last API call.
    """

    def __init__(self, calls_per_second: float):
        self.calls_per_second = calls_per_second
        self.last_call = 0

    async def wait(self):
        """
        Wait until the next API call can be made without exceeding the rate limit.
        """
        now = time.time()
        time_since_last_call = now - self.last_call
        if time_since_last_call < 1 / self.calls_per_second:
            await asyncio.sleep(1 / self.calls_per_second - time_since_last_call)
        self.last_call = time.time()


class ChemblAPIClient:
    """
       ChemblAPIClient handles interactions with the ChEMBL API, including fetching and retrying requests.

       Attributes:
           base_url (str): Base URL for the ChEMBL API.
           batch_size (int): Number of records to fetch per batch.
           timeout (int): Timeout for API requests.
           retries (int): Maximum number of retries for failed API requests.
           semaphore (asyncio.Semaphore): Semaphore to limit concurrent API requests.
           rate_limiter (RateLimiter): RateLimiter instance to manage API call rate.

    """

    def __init__(self, base_url: str):
        self.base_url = base_url
        self.batch_size = BATCH_SIZE
        self.timeout = TIMEOUT
        self.retries = MAX_RETRIES
        self.semaphore = asyncio.Semaphore(MAX_CONCURRENT_REQUESTS)
        self.rate_limiter = RateLimiter(API_RATE_LIMIT)

    @retry(stop=stop_after_attempt(MAX_RETRIES), wait=wait_exponential(multiplier=1, min=4, max=10))
    async def fetch_with_retry(self, session: aiohttp.ClientSession, url: str) -> Dict:
        """
        Fetch data from the given URL with retries on failure.

        Args:
            session (aiohttp.ClientSession): The aiohttp session used for making requests.
            url (str): The URL to fetch data from.

        Returns:
            Dict: Parsed JSON response from the API.

        Raises:
            aiohttp.ClientResponseError: If the response status indicates an error.
            json.JSONDecodeError: If the response cannot be parsed as JSON.
        """

        await self.rate_limiter.wait()
        async with self.semaphore:
            try:
                async with session.get(url, timeout=self.timeout) as response:
                    if response.status == 429:  # Too Many Requests
                        logger.warning('Rate limit hit. Waiting before retrying...')
                        await asyncio.sleep(60)  # Wait for 1 minute before retrying
                        raise aiohttp.ClientResponseError(response.request_info, response.history, status=429)
                    elif response.status == 500:
                        logger.error(f'Internal Server Error for URL {url}. Retrying...')
                        raise aiohttp.ClientResponseError(response.request_info, response.history, status=500)
                    response.raise_for_status()
                    data = await response.text()
                    try:
                        parsed_data = json.loads(data)
                        logger.debug(f'Successfully fetched and parsed data from URL {url}')
                        return parsed_data
                    except json.JSONDecodeError as json_err:
                        logger.error(f'Failed to parse JSON from response. URL: {url}. Error: {str(json_err)}')
                        logger.error(f'Raw data (first 200 chars): {data[:10]}')
                        raise
            except aiohttp.ClientResponseError as e:
                logger.error(f'HTTP error {e.status} fetching data from URL {url}: {e.message}')
                raise
            except aiohttp.ClientError as e:
                logger.error(f'Client error fetching data from URL {url}: {e}')
                raise
            except asyncio.TimeoutError:
                logger.error(f'Timeout error fetching data from URL {url}')
                raise

    async def fetch_total_records(self, session: aiohttp.ClientSession) -> int:
        """
         Fetch the total number of records available from the API.

         Args:
             session (aiohttp.ClientSession): The aiohttp session used for making requests.

         Returns:
             int: Total number of records available.

         Raises:
             Exception: If fetching total records fails.
         """
        url = f'{self.base_url}?limit=1'
        data = await self.fetch_with_retry(session, url)
        total_records = data['page_meta']['total_count']
        logger.info(f'Total records to fetch: {total_records}')
        return total_records

    async def fetch_data_chunk(self, session: aiohttp.ClientSession, offset: int = 0) -> Optional[Dict]:
        """
          Fetch a chunk of data from the API.

          Args:
              session (aiohttp.ClientSession): The aiohttp session used for making requests.
              offset (int): Offset for pagination.

          Returns:
              Optional[Dict]: Parsed JSON response from the API, or None if an error occurs.
        """
        url = f'{self.base_url}?limit={self.batch_size}&offset={offset}'
        try:
            data = await self.fetch_with_retry(session, url)
            logger.info(f'Fetched data chunk with offset {offset}')
            return data
        except aiohttp.ClientResponseError as e:
            if e.status == 429:
                logger.warning('Rate limit hit. Retrying after a delay...')
                await asyncio.sleep(60)  # Sleep to avoid immediate retry
            logger.error(f'HTTP error {e.status} fetching data with offset {offset}: {e.message}')
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            logger.error(f'Error fetching data with offset {offset}: {e}')
        return None

    async def fetch_data_chunks(self, session: aiohttp.ClientSession) -> List[List[Dict]]:
        """
        Fetch all data chunks from the API.

        Args:
            session (aiohttp.ClientSession): The aiohttp session used for making requests.

        Returns:
            List[List[Dict]]: List of data chunks, each containing a list of records.
        """
        total_records = await self.fetch_total_records(session)
        all_data = []
        tasks = []
        for offset in range(0, total_records, self.batch_size):
            tasks.append(asyncio.create_task(self.fetch_data_chunk(session, offset)))
        chunk_results = await asyncio.gather(*tasks)
        for chunk in chunk_results:
            if chunk is None:
                logger.error(f'Chunk is None')
                continue
            if chunk:
                records = chunk.get('chembl_id_lookups') or chunk.get('molecules', [])
                if isinstance(records, list):
                    all_data.append(records)
                else:
                    logger.error(f'Expected records to be a list, but got {type(records)}: ')
            else:
                logger.error(f'Expected chunk to be a dict, but got {type(chunk)}: ')
        return all_data
