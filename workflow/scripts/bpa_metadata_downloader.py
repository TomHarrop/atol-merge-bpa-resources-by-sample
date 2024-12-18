#!/usr/bin/env python3

import ckanapi
import os

# Globals
dp_url = "https://data.bioplatforms.com"


def get_apikey():
    apikey = os.getenv("BPI_APIKEY")
    if not apikey:
        raise ValueError("Set the BPI_APIKEY environment variable.")
    return apikey


def main():
    apikey = get_apikey()
    remote = ckanapi.RemoteCKAN(dp_url, apikey=apikey)
    all_packages = remote.action.package_list(
        rows=0, include_datasets=True, include_private=True, apikey=apikey
    )


main()

if __name__ == "__main__":

    main()
