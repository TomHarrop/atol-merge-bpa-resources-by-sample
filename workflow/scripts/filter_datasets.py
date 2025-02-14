#!/usr/bin/env/python3

from collections import Counter
from utils import setup_logging
import gzip
import json
import jsonlines
import csv


class ConfigSection:
    def __init__(self, field_names, accepted_values):
        self.field_names = field_names
        self.accepted_values = accepted_values


class Config:
    def __init__(self, config_file):
        with open(config_file, "r") as file:
            config_data = json.loads(file.read())
        for key, value in config_data.items():
            setattr(
                self, key, ConfigSection(value["field_names"], value["accepted_values"])
            )

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


def choose_value(package, keys_to_check, accepted_values):
    """
    Returns a tuple of (chosen_value, accept_value).

    Package is a dict parsed from json.

    keys_to_check is an ordered list.

    If package has any keys_to_check whose value is in accepted_values, the
    value of that keys_to_check is returned and accept_value is True.

    If the package has no keys_to_check whose value is in accepted_values, the
    value of the first keys_to_check is returned.

    If the package has keys_to_check at all, the value is None.

    """
    values = {key: get_nested_value(package, key) for key in keys_to_check}
    first_value = None
    first_key = None

    for key, value in values.items():
        if value is not None:
            if value in accepted_values:
                return (value, key, True)
            if first_value is None:
                first_value = value
                first_key = key

    return (first_value, first_key, False)


def get_nested_value(d, key):
    """
    Retrieve the value from a nested dictionary using a dot-notated key.
    """
    keys = key.split(".")
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return None
    return d


def read_json(json_file):
    """
    Each line is a Package.
    """
    with gzip.open(json_file, "rb") as f:
        data = jsonlines.Reader(f)
        for line in data:
            yield line


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
