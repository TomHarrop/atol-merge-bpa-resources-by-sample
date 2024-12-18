#!/usr/bin/env python3

import ckanapi
import os
import multiprocessing as mp

# Globals
dp_url = "https://data.bioplatforms.com"
threads = 8
num_records = 100  # in prod, set to 0


def get_apikey():
    apikey = os.getenv("BPI_APIKEY")
    if not apikey:
        raise ValueError("Set the BPI_APIKEY environment variable.")
    return apikey


def main():
    apikey = get_apikey()
    remote = ckanapi.RemoteCKAN(dp_url, apikey=apikey)
    all_packages = remote.action.package_list(
        rows=num_records, include_datasets=True, include_private=True, apikey=apikey
    )
    print(f"Found {len(all_packages)} packages.")
    raise ValueError()


main()

if __name__ == "__main__":

    main()
