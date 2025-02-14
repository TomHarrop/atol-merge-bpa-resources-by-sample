#!/usr/bin/env python3

from collections import Counter
from config_parser import Config
from metadata_parser import choose_value, get_nested_value
from utils import setup_logging, read_json, write_json


def main():
    config = Config(mapping_config_file)

    # track key usage
    key_usage_counters = {x: Counter() for x, y in config}
    parsed_key_counters = {x: Counter() for x, y in config}

    # this is one big dict of id: package
    data = next(read_json(json_file))

    id_set = set()

    i = 0

    for package in data.values():
        id = package["id"]

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
            parsed_key_counters[metadata_field].update([key])

    write_json(key_usage_counters, key_usage_counts_file)
    write_json(parsed_key_counters, parsed_key_counts_file)
    quit(1)


if __name__ == "__main__":

    # Globals
    try:
        logger = setup_logging(snakemake.log[0])

    except NameError as e:
        import tempfile

        logfile = tempfile.mkstemp()[1]
        print(f"Logging to {logfile}")

        logger = setup_logging(logfile)
        logger.warning("No Snakemake object. Running script as standalone.")

        json_file = "results/filtered_datasets.jsonl.gz"
        mapping_config_file = "config/organism_mapping_config.json"
        key_usage_counts_file = "test/key_usage_counts.jsonl.gz"
        parsed_key_counts_file = "test/parsed_key_counts.jsonl.gz"

    main()
