import logging

import pendulum
from airflow import DAG
from airflow.operators.empty import EmptyOperator
from airflow.operators.python import PythonOperator
from airflow.providers.postgres.operators.postgres import PostgresOperator
from botocore.exceptions import ClientError

from config import (
    S3_BUCKET_NAME,
    S3_INPUT_FILES,
    S3_FINGERPRINTS_FOLDER_PATH,
    S3_HOOK,
)
from populate_presentation_layer.populate_dim_table import POPULATE_DIM_TABLE
from populate_presentation_layer.populate_fact_table import \
    ingest_top_10_similarities
from processor.process_similarity_score import ProcessSimilarityScore
from processor.process_source_molecules import ProcessSourceMolecules
from processor.process_target_molecules import ProcessTargetMolecules
from utils.s3_operations import AsyncS3Saver
from utils.slack_notifier import slack_notification

logging.basicConfig(level=logging.INFO)


def check_for_new_files(ti, **kwargs):
    s3_hook = S3_HOOK
    execution_date = kwargs['execution_date']
    current_month = execution_date.strftime('%m')
    current_year = execution_date.strftime('%Y')
    prefix = f'{S3_INPUT_FILES}data_{current_month}_{current_year}'

    try:
        file_list = s3_hook.list_keys(
            bucket_name=S3_BUCKET_NAME,
            prefix=prefix,
        )
        if file_list:
            logging.info(f'Found {len(file_list)} files for {current_month}/{current_year}')
            ti.xcom_push(
                key='new_file_list',
                value=file_list,
            )
            return True
        else:
            logging.info(f'No files found for {current_month}/{current_year}')
            return False
    except ClientError as e:
        logging.error(f'Error checking for files: {e}')
        return False


def process_target_molecules(ti):
    file_list = ti.xcom_pull(
        task_ids='check_for_new_files',
        key='new_file_list',
    )
    if not file_list:
        logging.info("No new files to process.")
        return

    s3_saver = AsyncS3Saver()
    target_processor = ProcessTargetMolecules(
        s3_saver=s3_saver,
        input_folder=S3_INPUT_FILES,
    )
    target_mols_df = target_processor.load_and_clean_csv_files(file_list)

    ti.xcom_push(
        key='target_mols_df',
        value=target_mols_df,
    )


def process_source_molecules(ti):

    s3_saver = AsyncS3Saver()
    source_processor = ProcessSourceMolecules(
        s3_saver=s3_saver,
        s3_fingerprints_folder=S3_FINGERPRINTS_FOLDER_PATH,
    )
    fingerprints_df = source_processor.load_and_precompute_fingerprints()
    ti.xcom_push(
        key='fingerprints_df',
        value=fingerprints_df,
    )


def compute_similarity_scores(ti):
    s3_saver = AsyncS3Saver()
    target_mols_df = ti.xcom_pull(
        key='target_mols_df',
        task_ids='process_target_molecules',
    )
    fingerprints_df = ti.xcom_pull(
        key='fingerprints_df',
        task_ids='process_source_molecules',
    )
    process_similarity_score = ProcessSimilarityScore(
        s3_saver=s3_saver,
        target_mols_df=target_mols_df,
        source_mols_df=fingerprints_df,
    )
    top_10_similar_results = process_similarity_score.compute_similarity_scores()
    ti.xcom_push(
        key='top_10_similar_results',
        value=top_10_similar_results,
    )


def populate_fact_table(ti):
    top_10_similar_results = ti.xcom_pull(
        key='top_10_similar_results',
        task_ids='compute_similarity_scores',
    )

    for chembl_id, top_10_similar in top_10_similar_results:
        if not top_10_similar.empty:
            ingest_top_10_similarities(
                chembl_id,
                top_10_similar,
            )
        else:
            logging.warning(f'No top 10 similarities found for {chembl_id}')
    logging.info('All top 10 similarities ingested successfully')


with DAG(
        dag_id='chembl_molecule_processing',
        description='Monthly dag that computes similarity score between target and source molecules',
        start_date=pendulum.datetime(2024, 8, 1),
        schedule_interval='0 0 1 * *',
        tags=['Final_work'],
        catchup=False,
        on_failure_callback=slack_notification,
) as dag:
    start_op = EmptyOperator(
        task_id='start'
    )

    check_for_new_files_task = PythonOperator(
        task_id='check_for_new_files',
        python_callable=check_for_new_files,
    )

    process_target_molecules_task = PythonOperator(
        task_id='process_target_molecules',
        python_callable=process_target_molecules,
    )

    process_source_molecules_task = PythonOperator(
        task_id='process_source_molecules',
        python_callable=process_source_molecules,
    )

    compute_similarity_scores_task = PythonOperator(
        task_id='compute_similarity_scores',
        python_callable=compute_similarity_scores,
    )

    populate_fact_table_task = PythonOperator(
        task_id='populate_fact_table',
        python_callable=populate_fact_table,
    )

    populate_dim_table_task = PostgresOperator(
        task_id='populate_dim_molecules',
        postgres_conn_id='db_conn',
        sql=POPULATE_DIM_TABLE
    )

    end_op = EmptyOperator(
        task_id='end'
    )

    start_op >> check_for_new_files_task >> [
        process_target_molecules_task,
        process_source_molecules_task,
        end_op
    ]

    process_target_molecules_task >> compute_similarity_scores_task
    process_source_molecules_task >> compute_similarity_scores_task

    compute_similarity_scores_task >> populate_fact_table_task >> populate_dim_table_task >> end_op
