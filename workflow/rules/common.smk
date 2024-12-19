import os
from snakemake.logging import logger


def get_apikey():
    apikey = os.getenv("BPI_APIKEY")
    if not apikey:
        logger.warning("Set the BPI_APIKEY environment variable.")
        return None
    return apikey


configfile: "config/config.yaml"
