#!/usr/bin/env/python3

import jsonlines
import gzip

# Globals
json_file = "datasets.jsonl.gz"


def read_json(json_file):
    '''
    Each line is a Package.
    '''
    with gzip.open(json_file, "rb") as f:
        data = jsonlines.Reader(f)
        for line in data:
            yield line


data = read_json(json_file)

package  = next(data)


