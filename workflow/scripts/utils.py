#!/usr/bin/env python3

from snakemake.logging import logger
import csv
import gzip
import jsonlines
import logging


def get_jsonlines_output_handle(file):
    """
    Open a gzip file and return a jsonlines.Writer object.
    """
    f = gzip.open(file, "wb")
    return jsonlines.Writer(f)


def read_json(json_file):
    """
    Each line is a Package.
    """
    with gzip.open(json_file, "rb") as f:
        data = jsonlines.Reader(f)
        for line in data:
            yield line


def setup_logging(logfile, log_level="INFO"):
    log_level = log_level.upper()
    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)
    logger.logger.setLevel(log_level)
    return logger


def write_decision_log_to_csv(decision_log, file_path):
    """
    Write the decision log to a CSV file.
    """
    with open(file_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Write the header
        header = ["id"] + list(next(iter(decision_log.values())).keys())
        writer.writerow(header)
        # Write the rows
        for id, decisions in decision_log.items():
            row = [id] + list(decisions.values())
            writer.writerow(row)


def write_json(data, file):
    with gzip.open(file, "wb") as f:
        jsonlines.Writer(f).write(data)
