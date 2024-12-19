#!/usr/bin/env/python3

import jsonlines
import gzip


def setup_logging(logfile):
    import logging
    from snakemake.logging import logger

    file_handler = logging.FileHandler(logfile)
    logger.logfile_handler = file_handler
    logger.logger.addHandler(logger.logfile_handler)
    return logger


def read_json(json_file):
    """
    Each line is a Package.
    """
    with gzip.open(json_file, "rb") as f:
        data = jsonlines.Reader(f)
        for line in data:
            yield line


def write_json(data, file):
    with gzip.open(file, "wb") as f:
        jsonlines.Writer(f).write_all(data)


def main():
    data = read_json(json_file)

    sample_id_dict = {}
    bioplatforms_sample_id_dict = {}
    equivalent_samples = []

    for i, package in enumerate(data):
        bp_id_exists = False
        if "bioplatforms_sample_id" in package:
            bp_id_exists = True
            bp_id = package["bioplatforms_sample_id"]
            if bp_id in bioplatforms_sample_id_dict:
                bioplatforms_sample_id_dict[bp_id].append(package)
            else:
                bioplatforms_sample_id_dict[bp_id] = [package]
        if "sample_id" in package:
            my_id = package["sample_id"]
            if bp_id_exists:
                equivalent_samples.append((my_id, bp_id))
            if my_id in sample_id_dict:
                sample_id_dict[my_id].append(package)
            else:
                sample_id_dict[my_id] = [package]

    write_json(equivalent_samples, equivalent_samples_file)
    write_json(sample_id_dict, sample_id_file)
    write_json(bioplatforms_sample_id_dict, bioplatforms_sample_id_file)


if __name__ == "__main__":

    # Globals
    if snakemake:

        logfile = snakemake.log[0]
        logger = setup_logging(logfile)

        json_file = snakemake.input["json"]
        equivalent_samples_file = snakemake.output["equivalent_samples"]
        sample_id_file = snakemake.output["sample_id_dict"]
        bioplatforms_sample_id_file = snakemake.output["bioplatforms_id_dict"]

    else:

        json_file = "datasets.jsonl.gz"
        equivalent_samples_file = "equivalent_samples.jsonl.gz"
        sample_id_file = "sample_id.jsonl.gz"
        bioplatforms_sample_id_file = "bioplatforms_sample_id.jsonl.gz"

    main()
