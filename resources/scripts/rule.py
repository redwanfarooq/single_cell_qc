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


def get_path(wildcards, path: str | None) -> str | None:
    """
    Get path from unformatted path string by substituting sample name.

    Arguments:
        ``wildcards``: Snakemake ``wildcards`` object.
        ``path``: unformatted path string (use {sample} as placeholder for sample name).

    Returns:
        Formatted path string or ``None`` if path not provided.
    """

    return path.format(sample=wildcards.sample) if path is not None else None
