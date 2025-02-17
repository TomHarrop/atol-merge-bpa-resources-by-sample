#!/usr/bin/env/python3

from collections import Counter
from config_parser import Config
from metadata_parser import get_nested_value, choose_value
from utils import (
    get_jsonlines_output_handle,
    read_json,
    setup_logging,
    write_decision_log_to_csv,
    write_json,
)


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

    # output
    logger.info(f"Writing filtered datasets to {filtered_datasets_file}")
    output_writer = get_jsonlines_output_handle(filtered_datasets_file)

    # logging
    i = 0
    j = 0

    for package in data:
        id = package["id"]
        package_decisions = {}

        logger.debug(f"Processing package {id}")

        # There shouldn't be any duplicate IDs
        if id in id_set:
            raise ValueError(f"Duplicate ID: {id}")
        id_set.add(id)

        for metadata_field, config_section in config:
            logger.debug(f"Processing metadata field {metadata_field}")
            for field in config_section.field_names:
                if get_nested_value(package, field) is not None:
                    key_usage_counters[metadata_field].update([field])

            value, key, keep = choose_value(
                package, config_section.field_names, config_section.accepted_values
            )

            # Macrodata refiners hate this one weird key. This is a manual
            # override. If the package has no context_keys whose value is in
            # accepted_data_context, but it does have a key called
            # "genome_data" with value "yes", keep_package is True.
            if (
                metadata_field == "data_context"
                and "genome_data" in package
                and not keep
            ):
                logger.debug("Checking genome_data field")
                if package["genome_data"] == "yes":
                    logger.debug("Setting keep to True")
                    value, key, keep = ("yes", "genome_data", True)

            logger.debug(f"value: {value}, key: {key}, keep: {keep}")

            parsed_value_counters[metadata_field].update([value])
            parsed_key_counters[metadata_field].update([key])

            # record the decision for this metadata_field
            decision_key = f"{metadata_field}_accepted"
            package_decisions[decision_key] = keep
            package_decisions[metadata_field] = value
            logger.debug(f"Decision: {package_decisions}")

        # summarize the decisions for this package
        logger.debug(f"Summarizing decisions for package")
        package_decisions["package_accepted"] = all(package_decisions.values())
        decision_log[id] = package_decisions

        # keep the package if all metadata_fields are accepted
        if package_decisions["package_accepted"]:
            logger.debug(f"Keeping package {id}")
            output_writer.write({id: package})
            j += 1

        i += 1

    logger.info(f"Processed {i} packages")
    logger.info(f"Wrote {j} packages to {filtered_datasets_file}")

    logger.info("Writing count files")
    write_json(key_usage_counters, key_usage_counts_file)
    write_json(parsed_value_counters, parsed_value_counts_file)
    write_json(parsed_key_counters, parsed_key_counts_file)
    logger.info("Writing decision log")
    write_decision_log_to_csv(decision_log, decision_log_file)
    logger.info("Done")


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

        import argparse
        import tempfile

        def parse_args():
            parser = argparse.ArgumentParser(
                description="Filter datasets based on configuration."
            )
            input_group = parser.add_argument_group("Input files")
            output_group = parser.add_argument_group("Output files")
            options = parser.add_argument_group("Options")
            input_group.add_argument(
                "--json",
                type=str,
                default="resources/datasets.jsonl.gz",
                help="Path to the JSON file.",
                dest="json_file",
            )
            input_group.add_argument(
                "--filtering_config",
                type=str,
                default="config/dataset_filtering_config.json",
                help="Path to the filtering configuration file.",
                dest="filtering_config_file",
            )
            output_group.add_argument(
                "--filtered_datasets",
                type=str,
                default="test/filtered_datasets.jsonl.gz",
                help="Path to the output filtered datasets file.",
                dest="filtered_datasets_file",
            )
            output_group.add_argument(
                "--decision_log",
                type=str,
                default="test/decision_log.csv",
                help="Path to the decision log file.",
                dest="decision_log_file",
            )
            output_group.add_argument(
                "--key_usage_counts",
                type=str,
                default="test/key_usage_counts.jsonl.gz",
                help="Path to the key usage counts file.",
                dest="key_usage_counts_file",
            )
            output_group.add_argument(
                "--parsed_value_counts",
                type=str,
                default="test/parsed_value_counts.jsonl.gz",
                help="Path to the parsed value counts file.",
                dest="parsed_value_counts_file",
            )
            output_group.add_argument(
                "--parsed_key_counts",
                type=str,
                default="test/parsed_key_counts.jsonl.gz",
                help="Path to the parsed key counts file.",
                dest="parsed_key_counts_file",
            )
            options.add_argument(
                "--log_level",
                type=str,
                default="INFO",
                help="Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR, CRITICAL).",
                dest="log_level",
            )
            return parser.parse_args()

        logfile = tempfile.mkstemp()[1]
        print(f"Logging to {logfile}")

        args = parse_args()
        globals().update(vars(args))

        logger = setup_logging(logfile, log_level=log_level)
        logger.warning("No Snakemake object. Running script as standalone.")

    main()
