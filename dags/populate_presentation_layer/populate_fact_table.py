from typing import List, Tuple

import pandas as pd
from psycopg2 import sql
from psycopg2.extras import execute_batch

from config import POSTGRES_HOOK
from utils.logging_config import setup_logging

logger = setup_logging()


def insert_similarity_scores(data: List[Tuple[str, str, float, bool]]):
    """
    Insert similarity scores into the database.

    :param data: List of tuples (target_chembl_id, source_chembl_id, tanimoto_similarity, has_duplicates_of_last_largest_score)
    """
    insert_query = sql.SQL("""
        INSERT INTO fact_molecule_similarities 
        (target_chembl_id, source_chembl_id, tanimoto_similarity, has_duplicates_of_last_largest_score)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT (source_chembl_id, target_chembl_id) DO UPDATE
        SET tanimoto_similarity = EXCLUDED.tanimoto_similarity,
            has_duplicates_of_last_largest_score = EXCLUDED.has_duplicates_of_last_largest_score
    """)

    try:
        with POSTGRES_HOOK.get_conn() as conn:
            with conn.cursor() as cur:
                execute_batch(cur, insert_query, data)
            conn.commit()
    except Exception as e:
        logger.error(f'Error inserting data: {e}')
        raise


def ingest_top_10_similarities(target_chembl_id: str, top_10_similar: pd.DataFrame):
    """
    Ingest top 10 similar molecules for a target molecule into the database.

    :param target_chembl_id: ChEMBL ID of the target molecule
    :param top_10_similar: DataFrame containing top 10 similar molecules
    """
    try:
        data = [
            (
                target_chembl_id,
                row['chembl_id'],
                row['similarity_score'],
                row['has_duplicates_of_last_largest_score']
            )
            for _, row in top_10_similar.iterrows()
        ]
        insert_similarity_scores(data)
        logger.info(f'Successfully ingested top 10 similarities for target {target_chembl_id}')
    except Exception as e:
        logger.error(f'Error ingesting top 10 similarities for target {target_chembl_id}: {e}')
