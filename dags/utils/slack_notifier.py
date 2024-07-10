import json
from textwrap import dedent
from airflow.providers.slack.operators.slack_webhook import SlackWebhookOperator


def slack_notification(context):
    task_instance = context['task_instance']
    task = context['task']
    dag = context['dag']
    execution_date = context['execution_date']
    exception = context.get('exception')

    error_message = (
        f"Task {task.task_id} in DAG {dag.dag_id} failed on {execution_date}.\n"
        f"Exception: {str(exception)}\n"
        f"Task Parameters: {json.dumps(task_instance.params, indent=2)}\n"
        f"Logs: {task_instance.log_url}"
    )

    slack_msg = dedent(f"""
        :red_circle: *Task Failure Alert*
        *DAG*: `{dag.dag_id}`
        *Task*: `{task.task_id}`
        *Execution Time*: {execution_date}

        *Error Details*:
        ```
        {error_message}
        ```

        *Actions*:
        • Check the logs for more details: <{task_instance.log_url}|View Logs>
        • Manually rerun the task: <{task_instance.get_url()}|Rerun Task>
    """)

    failed_alert = SlackWebhookOperator(
        task_id='slack_notification',
        slack_webhook_conn_id='URL',
        message=slack_msg,
        channel="#project",
        username='Airflow Alert Bot',
        dag=dag
    )

    return failed_alert.execute(context=context)
