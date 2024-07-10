from multiprocessing import cpu_count

# etl.py constants
DATABASE_URL = 'URL'
CHEMBL_ID_LOOKUP_URL = 'https://www.ebi.ac.uk/chembl/api/data/chembl_id_lookup.json'
MOLECULE_URL = 'https://www.ebi.ac.uk/chembl/api/data/molecule.json'

# api_client constants
BATCH_SIZE = 1000
TIMEOUT = 300
MAX_RETRIES = 5
MAX_CONCURRENT_REQUESTS = 30
API_RATE_LIMIT = 50

# db_handler constants
API_CHUNK_SIZE = 100000

# fingerprint_computations constatnts
FINGERPRINT_CHUNK_SIZE = 150000
COMPUTE_NUM_PROCESSES = cpu_count
FPS_MOL_RADIUS = 2
FPS_BITS = 2048

# s3_operations constatnts
S3_BUCKET_NAME = 'S3_BUCKET_NAME'
S3_FINGERPRINTS_FOLDER_PATH = 'S3_FINGERPRINTS_FOLDER_PATH'

AWS_ACCESS_KEY = 'AWS_ACCESS_KEY'
AWS_SECRET_KEY = 'AWS_SECRET_KEY'
AWS_SESSION_TOKEN = 'AWS_SESSION_TOKEN'
