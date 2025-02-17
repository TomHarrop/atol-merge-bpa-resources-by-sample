#!/usr/bin/env python3

from collections import Counter
from config_parser import Config
from metadata_parser import choose_value, get_nested_value
from utils import setup_logging, read_json, write_json


def pick_values_for_package(package, config):
    section_metadata = {}
    for metadata_field, config_section in config:
        value, key, keep = choose_value(
            package, config_section.field_names, config_section.accepted_values
        )
        logger.debug(f"value: {value}, key: {key}, keep: {keep}")
        section_metadata[metadata_field] = value

    return section_metadata


def main():
    organism_config = Config(organism_mapping_config)
    sample_config = Config(sample_mapping_config)

    # track key usage
    key_usage_counters = {x: Counter() for x, y in organism_config}
    parsed_key_counters = {x: Counter() for x, y in organism_config}

    # this is one big dict of id: package
    data = read_json(json_file)

    id_set = set()

    for package_line in data:
        # check that we only got one item
        if len(package_line) != 1:
            raise ValueError(
                (f"Malformed jsonlines.\nExpected one item, got {len(package_line)}")
            )

        id, package = package_line.popitem()
        logger.debug(f"Processing package {id}")

        # There shouldn't be any duplicate IDs
        if id in id_set:
            raise ValueError(f"Duplicate ID: {id}")
        id_set.add(id)

        package_metadata = {
            "original_id": id,
            "organism": pick_values_for_package(package, organism_config),
            "sample": pick_values_for_package(package, sample_config),
        }

        print(package_metadata)
        quit(1)

        for metadata_field, config_section in organism_config:
            logger.debug(f"Processing metadata field {metadata_field}")
            for field in config_section.field_names:
                if get_nested_value(package, field) is not None:
                    key_usage_counters[metadata_field].update([field])

            value, key, keep = choose_value(
                package, config_section.field_names, config_section.accepted_values
            )

            logger.debug(f"value: {value}, key: {key}, keep: {keep}")
            parsed_key_counters[metadata_field].update([key])

    write_json(key_usage_counters, key_usage_counts_file)
    write_json(parsed_key_counters, parsed_key_counts_file)


if __name__ == "__main__":

    # Globals
    try:
        logger = setup_logging(snakemake.log[0])

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
                default="results/filtered_datasets.jsonl.gz",
                help="Path to the JSON file.",
                dest="json_file",
            )
            input_group.add_argument(
                "--organism_mapping_config",
                type=str,
                default="config/organism_mapping_config.json",
                help="Path to the filtering configuration file.",
                dest="organism_mapping_config",
            )
            input_group.add_argument(
                "--sample_mapping_config",
                type=str,
                default="config/sample_mapping_config.json",
                help="Path to the filtering configuration file.",
                dest="sample_mapping_config",
            )

            output_group.add_argument(
                "--key_usage_counts",
                type=str,
                default="test/key_usage_counts.jsonl.gz",
                help="Path to the key usage counts file.",
                dest="key_usage_counts_file",
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
