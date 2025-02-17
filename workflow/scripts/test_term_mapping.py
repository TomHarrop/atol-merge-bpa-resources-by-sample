#!/usr/bin/env python3

from config_parser import TermMapping

term_mapping_file = "config/term_mapping_list.json"

mappings = TermMapping(term_mapping_file)
print(mappings.data_context.accepted_values)

quit(1)
