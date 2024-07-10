import json
from multiprocessing import cpu_count

from airflow.hooks.base import BaseHook
from airflow.models import Variable
from airflow.providers.amazon.aws.hooks.s3 import S3Hook
from airflow.providers.postgres.hooks.postgres import PostgresHook

CHUNK_SIZE = int(Variable.get('CHUNK_SIZE'))
FPS_MOL_RADIUS = int(Variable.get('FPS_MOL_RADIUS'))
FPS_BITS = int(Variable.get('FPS_BITS'))

COMPUTE_NUM_PROCESSES = cpu_count()

POSTGRES_HOOK = PostgresHook(postgres_conn_id='db_conn')
S3_HOOK = S3Hook(aws_conn_id='aws_default')

S3_BUCKET_NAME = Variable.get('S3_BUCKET_NAME')
S3_FINGERPRINTS_FOLDER_PATH = Variable.get('S3_FINGERPRINTS_FOLDER_PATH')
S3_SIMILARITY_FOLDER_PATH = Variable.get('S3_SIMILARITY_FOLDER_PATH')
S3_TOP10_FOLDER_PATH = Variable.get('S3_TOP10_FOLDER_PATH')
S3_INPUT_FILES = Variable.get('S3_INPUT_FILES')

conn = BaseHook.get_connection('aws_default')
extra_config = json.loads(conn.extra or '{}')

AWS_ACCESS_KEY = conn.login
AWS_SECRET_KEY = conn.password
AWS_SESSION_TOKEN = extra_config.get('aws_session_token')
