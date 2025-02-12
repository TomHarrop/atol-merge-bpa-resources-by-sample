#!/usr/bin/env/python3

from collections import Counter
from utils import setup_logging
import gzip
import json
import jsonlines


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
        jsonlines.Writer(f).write(data)


def pick_id(package):
    for field in sample_name_fields:
        if field in package:
            return package[field]


def get_bioplatforms_project(package):
    try:
        return package["bioplatforms_project"]
    except KeyError:
        return None


def get_data_context(package, context_keys, accepted_data_contexts):
    my_context = None
    for key in context_keys:
        if key in package:
            my_context = package[key]
            break

    if my_context in accepted_data_contexts:
        return "atol_assembly"
    elif "genome_data" in package:
        if package["genome_data"] == "yes":
            return "atol_assembly"
    else:
        return my_context


def read_data_mapping_config(data_mapping_config_file):
    with open(data_mapping_config_file, "rb") as file:
        return json.loads(file.read())


def main():
    data_mapping_config = read_data_mapping_config(data_mapping_config_file)
    globals().update(data_mapping_config)

    logger.info(f"Parsing JSON file {json_file}")
    data = read_json(json_file)

    counters = {
        "bioplatforms_project": Counter(),
        "data_context": Counter(),
        "project_context": Counter(),
    }

    for package in data:
        id = pick_id(package)

        # TODO: count potential record identifiers

        try:
            package["dataset_id"]
        except KeyError as e:
            print(id)
            raise e


        bioplatforms_project = get_bioplatforms_project(package)
        counters["bioplatforms_project"].update([bioplatforms_project])

        data_context = get_data_context(package, context_keys, accepted_data_context)
        counters["data_context"].update([data_context])

        counters["project_context"].update([f"{bioplatforms_project}_{data_context}"])

        if (
            bioplatforms_project in accepted_bioplatforms_project
            and data_context == "atol_assembly"
        ):
            print(id)
            quit(1)

    print(counters)
    quit(1)


if __name__ == "__main__":

    # Globals
    try:
        logger = setup_logging(snakemake.log[0])

        json_file = snakemake.input["json"]
        data_mapping_config_file = snakemake.input["data_mapping_config"]
        rejected_samples_file = snakemake.output["rejected_samples"]

    except NameError as e:
        import tempfile

        logfile = tempfile.mkstemp()[1]
        print(f"Logging to {logfile}")

        logger = setup_logging(logfile)
        logger.warning("No Snakemake object. Running script as standalone.")

        json_file = "resources/datasets.jsonl.gz"
        data_mapping_config_file = "config/data_mapping_config.json"
        rejected_samples_file = "test/rejected_samples.csv"

    main()
