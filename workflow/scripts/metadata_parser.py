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
    values = {key: get_nested_value(package, key) for key in keys_to_check}
    first_value = None
    first_key = None

    for key, value in values.items():
        if value is not None:
            if value in accepted_values:
                return (value, key, True)
            if first_value is None:
                first_value = value
                first_key = key

    return (first_value, first_key, False)



def get_nested_value(d, key):
    """
    Retrieve the value from a nested dictionary using a dot-notated key.
    """
    keys = key.split(".")
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return None
    return d
