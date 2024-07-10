import logging


def setup_logging():
    logging.basicConfig(
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S',
        level=logging.INFO,
    )
    return logging.getLogger(__name__)
