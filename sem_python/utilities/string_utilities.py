from .error_utilities import error


def stripStringFromRight(original_string, string_to_strip):
    if original_string.endswith(string_to_strip):
        return original_string[:-len(string_to_strip)]
    else:
        return original_string
