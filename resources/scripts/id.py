"""
Functions for extracting/generating unique IDs for use with single cell data QC pipeline scripts.
"""


def paste(*args, sep: str = "") -> list[str]:
    """
    Vectorised string concatenation.

    Arguments:
        ``*args``: 2 or more lists of strings or data type coercible to string.\n
        ``sep``: String used as separator during concatenation. Default: ``""``.

    Returns:
        List of concetenated strings.
    """
    combs = zip(*args)
    return [sep.join(str(j) for j in i) for i in combs]


def sample_id(donor_id: list[str], hash_id: list[str]) -> list[str]:
    """
    Generate unique sample IDs from donor ID and hash ID.

    Arguments:
        ``donor_id``: List of strings specifying donor IDs.\n
        ``hash_id``: List of strings specifying hash IDs.

    Returns:
        List of unique sample IDs.
    """
    return paste(donor_id, hash_id, sep="-")
