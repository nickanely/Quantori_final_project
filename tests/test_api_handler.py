import time
from unittest.mock import patch, AsyncMock

import aiohttp
import pytest
from script.build_dwh.utils.config import CHEMBL_ID_LOOKUP_URL, API_RATE_LIMIT
from script.build_dwh.services.api_client import ChemblAPIClient, RateLimiter


@pytest.mark.asyncio
async def test_rate_limiter():
    rate_limiter = RateLimiter(1)
    start_time = time.time()
    await rate_limiter.wait()
    assert time.time() - start_time < 1


@pytest.mark.asyncio
async def test_fetch_with_retry():
    client = ChemblAPIClient(CHEMBL_ID_LOOKUP_URL, API_RATE_LIMIT)
    with patch('aiohttp.ClientSession.get', new_callable=AsyncMock) as mock_get:
        mock_get.return_value.__aenter__.return_value.status = 200
        mock_get.return_value.__aenter__.return_value.text = AsyncMock(
            return_value='{"page_meta": {"total_count": 100}}')
        async with aiohttp.ClientSession() as session:
            response = await client.fetch_with_retry(session, 'http://test-url')
            assert response['page_meta']['total_count'] == 100


@pytest.mark.asyncio
async def test_fetch_total_records():
    client = ChemblAPIClient(CHEMBL_ID_LOOKUP_URL, API_RATE_LIMIT)
    with patch.object(client, 'fetch_with_retry', new_callable=AsyncMock) as mock_fetch:
        mock_fetch.return_value = {'page_meta': {'total_count': 100}}
        async with aiohttp.ClientSession() as session:
            total_records = await client.fetch_total_records(session)
            assert total_records == 100
