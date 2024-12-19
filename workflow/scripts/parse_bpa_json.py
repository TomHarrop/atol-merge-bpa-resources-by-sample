#!/usr/bin/env/python3

import jsonlines
import gzip

# Globals
json_file = "datasets.jsonl.gz"


def read_json(json_file):
    """
    Each line is a Package.
    """
    with gzip.open(json_file, "rb") as f:
        data = jsonlines.Reader(f)
        for line in data:
            yield line


data = read_json(json_file)

sample_id_dict = {}
bioplatforms_sample_id_dict = {}
equivalent_samples = []

for i, package in enumerate(data):
    bp_id = False
    if "bioplatforms_sample_id" in package:
        bp_id = True
        my_id = package["bioplatforms_sample_id"]
        if my_id in bioplatforms_sample_id_dict:
            bioplatforms_sample_id_dict[my_id].append(package)
        else:
            bioplatforms_sample_id_dict[my_id] = [package]
    if "sample_id" in package:
        my_id = package["sample_id"]
        if bp_id:
            equivalent_samples.append((sample_id, bioplatforms_sample_id))
        if my_id in sample_id_dict:
            sample_id_dict[my_id].append(package)
        else:
            sample_id_dict[my_id] = [package]

equivalent_samples
