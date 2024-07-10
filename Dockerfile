# Dockerfile
FROM apache/airflow:2.8.1

# Set environment variables
ENV AIRFLOW_HOME=/opt/airflow

# Copy the requirements file and install dependencies
COPY dags/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the DAGs directory
COPY dags /opt/airflow/dags

# Set the user to airflow
USER airflow
