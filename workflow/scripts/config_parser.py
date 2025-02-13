#!/usr/bin/env python3

import json


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
