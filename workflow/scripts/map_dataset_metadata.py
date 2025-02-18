#!/usr/bin/env python3

from collections import Counter
from config_parser import Config, TermMapping
from metadata_parser import choose_value, get_nested_value
from utils import setup_logging, read_json, write_json, get_jsonlines_output_handle


def pick_values_for_package(package, config, term_mapping, log_writer, counters=None):
    section_metadata = {}
    for metadata_field, config_section in config:
        # count them if we can
        if counters is not None:
            for field in config_section.field_names:
                if get_nested_value(package, field) is not None:
                    counters["field_present"][metadata_field].update([field])

        value, key, keep = choose_value(
            package, config_section.field_names, config_section.accepted_values
        )
        logger.debug(f"value: {value}, key: {key}, keep: {keep}")

        mapped_value, log = map_term(metadata_field, value, term_mapping)
        log_writer.write({package["id"]: log})
        section_metadata[metadata_field] = mapped_value

        if counters is not None:
            counters["key_usage"][metadata_field].update([key])
            counters["value_usage"][metadata_field].update([value])

    return section_metadata


def get_counters(config):
    return {
        "field_present": {x: Counter() for x, y in config},
        "key_usage": {x: Counter() for x, y in config},
        "value_usage": {x: Counter() for x, y in config},
    }


def map_term(field_name, value, term_mapping):
    try:
        mapped_value = term_mapping.lookup(field_name, value)
    except KeyError:
        logger.debug(f"No mapping for {field_name}: {value}")
        mapped_value = value
    log = (field_name, value, mapped_value)
    return mapped_value, log


def main():
    term_mapping = TermMapping(term_mapping_config)

    organism_config = Config(organism_mapping_config)
    sample_config = Config(sample_mapping_config)
    reads_config = Config(reads_mapping_config)

    # track key usage
    organism_counters = get_counters(organism_config)
    sample_counters = get_counters(sample_config)
    reads_counters = get_counters(reads_config)

    # this is one big dict of id: package
    logger.info(f"Parsing JSON file {json_file}")
    data = read_json(json_file)

    logger.info(f"Writing modified metadata to {mapped_metadata_file}")
    output_writer = get_jsonlines_output_handle(mapped_metadata_file)
    term_mapping_log_writer = get_jsonlines_output_handle(term_mapping_log_file)

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
            "organism": pick_values_for_package(
                package,
                organism_config,
                term_mapping,
                term_mapping_log_writer,
                organism_counters,
            ),
            "sample": pick_values_for_package(
                package,
                sample_config,
                term_mapping,
                term_mapping_log_writer,
                sample_counters,
            ),
            "reads": pick_values_for_package(
                package,
                reads_config,
                term_mapping,
                term_mapping_log_writer,
                reads_counters,
            ),
        }

        output_writer.write(package_metadata)

    write_json(organism_counters, organism_counts_file)
    write_json(sample_counters, sample_counts_file)
    write_json(reads_counters, read_counts_file)


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

            # inputs

            input_group = parser.add_argument_group("Input files")

            input_group.add_argument(
                "--json",
                type=str,
                default="results/filtered_datasets.jsonl.gz",
                help="Path to the JSON file.",
                dest="json_file",
            )

            input_group.add_argument(
                "--term_mapping_config",
                type=str,
                default="config/term_mapping_list.json",
                help="Path to the filtering configuration file.",
                dest="term_mapping_config",
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

            input_group.add_argument(
                "--reads_mapping_config",
                type=str,
                default="config/reads_mapping_config.json",
                help="Path to the filtering configuration file.",
                dest="reads_mapping_config",
            )

            # outputs
            output_group = parser.add_argument_group("Output files")
            output_group.add_argument(
                "--mapped_metadata",
                type=str,
                default="test/mapped_metadata.jsonl.gz",
                help="Path to the key usage counts file.",
                dest="mapped_metadata_file",
            )
            output_group.add_argument(
                "--organism_counts",
                type=str,
                default="test/organism_counts.jsonl.gz",
                help="Path to the key usage counts file.",
                dest="organism_counts_file",
            )
            output_group.add_argument(
                "--sample_counts",
                type=str,
                default="test/sample_counts.jsonl.gz",
                help="Path to the parsed key counts file.",
                dest="sample_counts_file",
            )
            output_group.add_argument(
                "--read_counts",
                type=str,
                default="test/read_counts.jsonl.gz",
                help="Path to the parsed key counts file.",
                dest="read_counts_file",
            )
            output_group.add_argument(
                "--term_mapping_log",
                type=str,
                default="test/term_mapping_log.jsonl.gz",
                help="Path to the parsed key counts file.",
                dest="term_mapping_log_file",
            )

            # options
            options = parser.add_argument_group("Options")
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
