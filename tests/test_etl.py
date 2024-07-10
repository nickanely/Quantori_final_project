from unittest.mock import patch, AsyncMock

import pytest

from dags.script.build_dwh.etl import process_and_insert_data


@pytest.mark.asyncio
async def test_process_and_insert_data():
    mock_api_client = AsyncMock()
    mock_db_handler = AsyncMock()
    with patch('aiohttp.ClientSession', new_callable=AsyncMock):
        await process_and_insert_data(mock_api_client, AsyncMock())
        assert mock_db_handler.process_data.call_count > 0
