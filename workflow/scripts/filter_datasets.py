#!/usr/bin/env/python3

from collections import Counter
from utils import setup_logging
import gzip
import json
import jsonlines


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

    def get_csv_header(self):
        return ",".join(f"{x},keep_{x}" for x, y in self)


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
    my_value = None
    first_value = None
    for key in keys_to_check:
        if key in package:
            my_value = package[key]
            if my_value in accepted_values:
                return (my_value, True)
            if not first_value:
                first_value = my_value

    return (first_value, False)


def main():
    config = Config(filtering_config_file)

    logger.info(f"Parsing JSON file {json_file}")
    data = read_json(json_file)

    # track key usage
    key_usage_counters = {x: Counter() for x, y in config}
    value_counters = {x: Counter() for x, y in config}

    packages_to_keep = {}
    id_set = set()
    decision_log = [f"id,{config.get_csv_header()}"]

    for package in data:
        keep_dataset = False
        id = package["id"]

        # There shouldn't be any duplicate IDs
        if id in id_set:
            raise ValueError(f"Duplicate ID: {id}")
        id_set.add(id)

        # check which of the allowed_fields have been used in this package
        for metadata_field, config_section in config:
            for field in config_section.field_names:
                if field in package:
                    key_usage_counters[metadata_field].update([field])

            value, keep = choose_value(
                package, config_section.field_names, config_section.accepted_values
            )

            # Need a manual override for this one weird key. If the package has no
            # context_keys whose value is in accepted_data_context, but it does
            # have a key called "genome_data" with value "yes", keep_package is
            # True.
            if (
                metadata_field == "data_context"
                and "genome_data" in package
                and not keep
            ):
                if package["genome_data"] == "yes":
                    keep_context = True

            print(value)
            print(keep)
            quit(1)

        counters["data_context"].update([data_context])

        platform, keep_platform = choose_value(
            package, platform_keys, accepted_platform
        )
        counters["platform"].update([platform])

        if all([keep_organization, keep_context, keep_platform]):
            packages_to_keep[id] = package
            keep_dataset = True

        decision_log.append(
            f"{id},{keep_dataset},{organization_name},{keep_organization},{data_context},{keep_context},{platform},{keep_platform}"
        )

    write_json(packages_to_keep, filtered_datasets_file)
    write_json(counters, counters_file)
    with open(decision_log_file, "w") as f:
        f.write("\n".join(decision_log))


if __name__ == "__main__":

    # Globals
    try:
        logger = setup_logging(snakemake.log[0])

        json_file = snakemake.input["json"]
        filtering_config_file = snakemake.input["data_mapping_config"]
        filtered_datasets_file = snakemake.output["filtered_datasets"]
        decision_log_file = snakemake.output["decision_log"]
        counters_file = snakemake.output["counters"]

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
        counters_file = "test/key_counts.jsonl.gz"

    main()
