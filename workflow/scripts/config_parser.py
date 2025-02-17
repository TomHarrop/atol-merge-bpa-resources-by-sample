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
            field_names = value["field_names"]
            accepted_values = (
                value["accepted_values"] if "accepted_values" in value else None
            )
            setattr(self, key, ConfigSection(field_names, accepted_values))

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


class TermMapping(dict):

    def __init__(self, config_file):
        super().__init__()
        with open(config_file, "r") as file:
            config_data = json.loads(file.read())

        for field_name, mapping in config_data.items():
            self[field_name] = {}
            for output_value, accepted_values in mapping.items():
                for value in accepted_values:
                    self[field_name][value] = output_value


    def lookup(self, field_name, value):
        return self[field_name][value]
