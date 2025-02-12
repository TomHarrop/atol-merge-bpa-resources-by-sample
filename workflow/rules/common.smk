#!/usr/bin/env python3

import os
from functools import cache
from snakemake.logging import logger


@cache
def get_container(container_name):
    containers = config["containers"]
    if container_name not in containers:
        raise ValueError(f"Container {container_name} not found in config.")
    my_container = containers[container_name]
    return f"{my_container["prefix"]}://{my_container["url"]}:{my_container["tag"]}"


def get_apikey():
    apikey = os.getenv("BPI_APIKEY")
    if not apikey:
        logger.error("Set the BPI_APIKEY environment variable.")
        return None
    return apikey


configfile: "config/config.yaml"
globals().update(config)

