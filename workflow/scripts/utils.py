#!/usr/bin/env python3

import logging
from snakemake.logging import logger


def setup_logging(logfile):

    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)
    return logger
