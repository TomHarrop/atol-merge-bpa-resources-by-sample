#!/usr/bin/env/python3

from collections import Counter
from config_parser import Config
from metadata_parser import get_nested_value, choose_value
from utils import setup_logging, read_json, write_json, write_decision_log_to_csv


def main():
    config = Config(filtering_config_file)

    logger.info(f"Parsing JSON file {json_file}")
    data = read_json(json_file)

    # track key usage
    key_usage_counters = {x: Counter() for x, y in config}
    parsed_value_counters = {x: Counter() for x, y in config}
    parsed_key_counters = {x: Counter() for x, y in config}

    decision_log = {}
    packages_to_keep = {}
    id_set = set()

    for package in data:
        id = package["id"]
        package_decisions = {}

        # There shouldn't be any duplicate IDs
        if id in id_set:
            raise ValueError(f"Duplicate ID: {id}")
        id_set.add(id)

        for metadata_field, config_section in config:
            for field in config_section.field_names:
                if get_nested_value(package, field) is not None:
                    key_usage_counters[metadata_field].update([field])

            value, key, keep = choose_value(
                package, config_section.field_names, config_section.accepted_values
            )

            # Need a manual override for this one weird key. If the package has
            # no context_keys whose value is in accepted_data_context, but it
            # does have a key called "genome_data" with value "yes",
            # keep_package is True.
            if (
                metadata_field == "data_context"
                and "genome_data" in package
                and not keep
            ):
                if package["genome_data"] == "yes":
                    value, key, keep = ("yes", "genome_data", True)

            parsed_value_counters[metadata_field].update([value])
            parsed_key_counters[metadata_field].update([key])

            # record the decision for this metadata_field
            decision_key = f"{metadata_field}_accepted"
            package_decisions[decision_key] = keep
            package_decisions[metadata_field] = value

        # summarize the decisions for this package
        package_decisions["package_accepted"] = all(package_decisions.values())
        decision_log[id] = package_decisions

        # keep the package if all metadata_fields are accepted
        if package_decisions["package_accepted"]:
            packages_to_keep[id] = package

    write_json(packages_to_keep, filtered_datasets_file)
    write_json(key_usage_counters, key_usage_counts_file)
    write_json(parsed_value_counters, parsed_value_counts_file)
    write_json(parsed_key_counters, parsed_key_counts_file)
    write_decision_log_to_csv(decision_log, decision_log_file)


if __name__ == "__main__":

    # Globals
    try:
        logger = setup_logging(snakemake.log[0])

        json_file = snakemake.input["json"]
        filtering_config_file = snakemake.input["data_mapping_config"]
        filtered_datasets_file = snakemake.output["filtered_datasets"]
        decision_log_file = snakemake.output["decision_log"]
        key_usage_counts_file = snakemake.output["key_usage_counts"]
        parsed_value_counts_file = snakemake.output["parsed_value_counts"]
        parsed_key_counts_file = snakemake.output["parsed_key_counts"]

    except NameError as e:
        import tempfile

        logfile = tempfile.mkstemp()[1]
        print(f"Logging to {logfile}")

        logger = setup_logging(logfile)
        logger.warning("No Snakemake object. Running script as standalone.")

        json_file = "resources/datasets.jsonl.gz"
        filtering_config_file = "config/dataset_filtering_config.json"
        filtered_datsets_file = "test/filtered_datasets.jsonl.gz"
        decision_log_file = "test/decision_log.csv"
        key_usage_counts_file = "test/key_usage_counts.jsonl.gz"
        parsed_value_counts_file = "test/parsed_value_counts.jsonl.gz"
        parsed_key_counts_file = "test/parsed_key_counts.jsonl.gz"

    main()
