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


# TODO: this has been generalised to recurse through key sets, generalise the
# name and text.
def get_data_context(package, context_keys, accepted_data_contexts):
    """
    Returns a tuple of (context_value, keep_package).

    context_keys is an ordered list.

    If package has any context_keys whose value is in accepted_data_contexts,
    the value of that context_key is returned and keep_package is True.

    If the package has no context_keys whose value is in
    accepted_data_contexts, but it does have a key called "genome_data" with
    value "yes", keep_package is True.

    If the package has no context_keys whose value is in
    accepted_data_contexts, the value of the first context_key is returned.

    If the package has context_keys at all, the value is  None.

    """
    my_context = None
    first_context = None
    for key in context_keys:
        if key in package:
            my_context = package[key]
            if my_context in accepted_data_contexts:
                return (my_context, True)
            if not first_context:
                first_context = my_context

    return (first_context, False)


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
        "platform_key_usage": Counter(),
        "platform": Counter(),
        "project_and_context": Counter(),
        "package_identifier_usage": Counter(),
    }

    packages_to_keep = {}
    id_set = set()
    decision_log = [
        "id,bioplatforms_project,keep_project,data_context,keep_context,platform,keep_platform"
    ]

    for package in data:
        id = package["id"]

        # There shouldn't be any duplicate IDs
        if id in id_set:
            raise ValueError(f"Duplicate ID: {id}")
        id_set.add(id)

        for field in context_keys:
            if field in package:
                counters["context_key_usage"].update([field])

        for field in platform_keys:
            if field in package:
                counters["platform_key_usage"].update([field])

        bioplatforms_project, keep_project = get_data_context(
            package, ["bioplatforms_project"], accepted_bioplatforms_project
        )
        counters["bioplatforms_project"].update([bioplatforms_project])

        data_context, keep_context = get_data_context(
            package, context_keys, accepted_data_context
        )
        # need a manual override for this one weird key
        if "genome_data" in package:
            if package["genome_data"] == "yes":
                keep_context = True

        counters["data_context"].update([data_context])

        counters["project_and_context"].update(
            [f"{bioplatforms_project}_{data_context}"]
        )

        # hack for now - it might work with the existing get_data_context
        # function
        platform, keep_platform = get_data_context(
            package, platform_keys, accepted_platform
        )
        counters["platform"].update([platform])

        decision_log.append(
            f"{id},{bioplatforms_project},{keep_project},{data_context},{keep_context},{platform},{keep_platform}"
        )
        if keep_project and keep_context and keep_platform:
            packages_to_keep[id] = package

    write_json(packages_to_keep, filtered_datasets_file)
    with open(decision_log_file, "w") as f:
        f.write("\n".join(decision_log))


if __name__ == "__main__":

    # Globals
    try:
        logger = setup_logging(snakemake.log[0])

        json_file = snakemake.input["json"]
        data_mapping_config_file = snakemake.input["data_mapping_config"]
        filtered_datasets_file = snakemake.output["filtered_datasets"]
        decision_log_file = snakemake.output["decision_log"]

    except NameError as e:
        import tempfile

        logfile = tempfile.mkstemp()[1]
        print(f"Logging to {logfile}")

        logger = setup_logging(logfile)
        logger.warning("No Snakemake object. Running script as standalone.")

        json_file = "resources/datasets.jsonl.gz"
        data_mapping_config_file = "config/data_mapping_config.json"
        filtered_datsets_file = "test/filtered_datasets.jsonl.gz"
        decision_log_file = "test/decision_log.csv"

    main()
