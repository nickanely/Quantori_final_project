# ChEMBL ETL Pipeline

### Table of Contents
* Overview
* Dependencies
* Setup and Installation
* ETL Process
* Similarity Computation and Data Mart
* Data Warehouse Architecture
* Database views

## Overview
This project involves a complex ETL pipeline designed to fetch data from the ChEMBL database, process it, and compute molecular fingerprints and similarity scores. The processed data is then stored in an S3 bucket and ingested into a data mart, which is managed using Apache Airflow.

## Dependencies
* Python 3.8+
* Apache Airflow
* aiohttp
* aioboto3
* pandas
* rdkit
* s3fs
* pyarrow
* boto3
* tenacity

## Setup and Instalation

1. clone repository:
`git clone <repository_url>
cd <repository_directory>`
2. Create a virtual environment and activate it:
`python3 -m venv venv
source venv/bin/activate`
3. Install the required dependencies:
`pip install -r requirements.txt`

## ETL Process
The ETL process involves the following steps:

### Fetch Data from ChEMBL API:

* etl.py orchestrates the process of fetching data from the ChEMBL API using the ChemblAPIClient class.
* Data is processed and inserted into PostgreSQL using the DatabaseHandler class.
### Compute Fingerprints:
* Fingerprints for molecules are computed using RDKit and stored in the staging_compound_structures table.


## Similarity Computation and Data Mart

### Read Files from S3:
* Compound structures are read from S3 and used to compute Tanimoto similarity scores with other ChEMBL molecules.

### Compute Similarity Scores:

* The similarity scores for each molecule are computed using ProcessSimilarityScore.

### Save Results:
* Full similarity score tables are saved as parquet files and uploaded to S3.
* The top 10 most similar molecules are identified, and a flag is set for duplicates.

### Data Mart

* Dimension and fact tables are designed to store molecule properties and similarity scores.



## Data Warehouse Architecture

The data warehouse architecture is designed to support efficient storage, retrieval, and analysis of ChEMBL molecular data. The architecture is divided into three main layers: Staging, Data Mart, and Presentation.

### Layers

### 1. Staging Layer:

* Raw data from ChEMBL API is fetched and stored in PostgreSQL staging tables.
* The data includes compound structures and other related molecular information.

### 2. Data Mart Layer:

* The data mart is designed following a star schema to optimize query performance.
* Dimension Table:
`dim_molecules`: Contains details about each molecule, including properties and identifiers.
* Fact Table:
`fact_molecule_similairties`: Stores similarity scores between molecules along with metadata such as computation date and source.

### 3. Presentation Layer:

* This layer is designed to support end-user queries and reporting.
* It includes aggregated and processed data from the data mart, ready for business intelligence tools and dashboards.



### Database views:

* Average similarity score per source molecule :

![7-A.png](assets%2F7-A.png)

* Average deviation of alogp :

![7-B.png](assets%2F7-B.png)

* 10 random molecules similairty:

![8-A.png](assets%2F8-A.png)

* Next most similar molecule:
![8-B.png](assets%2F8-B.png)

* Average grouping:
![8-C.png](assets%2F8-C.png)
![8-C-1.png](assets%2F8-C-1.png)