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
    for field in package_identifier_fields:
        if field in package:
            return package[field]


def get_bioplatforms_project(package):
    try:
        return package["bioplatforms_project"]
    except KeyError:
        return None


def get_data_context(
    package, context_keys, accepted_data_contexts, assembly_value="atol_assembly"
):
    """
    context_keys is an ordered list. If package has any context_keys whose
    value is in accepted_data_contexts, the function returns assembly_value.
    Otherwise, it returns the value of the first context_key that was present
    in the package. If no context_keys are present, it returns None.
    """
    my_context = None
    first_context = None
    for key in context_keys:
        if key in package:
            my_context = package[key]
            if my_context in accepted_data_contexts:
                return assembly_value
            if not first_context:
                first_context = my_context

    if "genome_data" in package:
        if package["genome_data"] == "yes":
            return assembly_value
    else:
        return first_context


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
        "context_key_usage": Counter(),
        "data_context": Counter(),
        "project_and_context": Counter(),
        "package_identifier_usage": Counter(),
    }

    for package in data:
        id = package["id"]

        for field in context_keys:
            if field in package:
                counters["context_key_usage"].update([field])

        bioplatforms_project = get_bioplatforms_project(package)
        counters["bioplatforms_project"].update([bioplatforms_project])

        data_context = get_data_context(package, context_keys, accepted_data_context)
        counters["data_context"].update([data_context])

        counters["project_and_context"].update(
            [f"{bioplatforms_project}_{data_context}"]
        )

        if bioplatforms_project in accepted_bioplatforms_project and not data_context:
            print(id)
            quit(1)

    print(counters["data_context"])
    print(counters["context_key_usage"])
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
