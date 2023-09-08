"""
Functions for use in Snakemake rule definitions.
"""

import yaml


def parse_info(info: dict) -> dict:
    """
    Parse dictionary of sample info.

    Arguments:
        ``info``: dictionary of sample info.

    Returns:
        Dictionary of parsed sample info.
    """
    samples = list(info.keys())

    return {"samples": samples}


def get_hto_metadata(wildcards, info: dict) -> str:
    """
    Get cell hashing metadata string.

    Arguments:
        ``wildcards``: Snakemake ``wildcards`` object.
        ``info``: dictionary of sample info.

    Returns:
        YAML-formatted cell hashing metadata string.
    """
    return yaml.dump(info[wildcards.sample])
