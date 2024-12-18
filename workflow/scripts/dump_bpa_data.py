#!/usr/bin/env python3

import ckanapi
import os

# ckanapi search datasets -a $BPA_APIKEY include_private=true -O datasets.jsonl.gz -z  -r https://data.bioplatforms.com
# Globals
dp_url = "https://data.bioplatforms.com"


def get_apikey():
    apikey = os.getenv("BPI_APIKEY")
    if not apikey:
        raise ValueError("Set the BPI_APIKEY environment variable.")
    return apikey


def main():
    pass


apikey = get_apikey()
remote = ckanapi.RemoteCKAN(dp_url, apikey=apikey)
all_results = remote.action.package_search(include_private=True)

x = remote.call_action("current_package_list")

len(x["results"])

all_results = remote.call_action(

    "package_search", {"include_private": True, "rows": 10000}
)

print(f"Found {len(all_packages)} packages.")
raise ValueError()
