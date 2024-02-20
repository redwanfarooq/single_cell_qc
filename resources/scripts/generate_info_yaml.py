#!/bin/env python


"""
Generates info YAML for use with single cell data QC pipeline.
Requires:
- Metadata CSV file with the following fields:
    donor: donor ID
    pool: pool ID
    sample: sample name
    hto: hashtag antibody ID (or 'None' if not used)
    cells_loaded: number of cells loaded
"""


# ==============================
# MODULES
# ==============================
import os
import yaml
import docopt
from loguru import logger
import pandas as pd
from id import sample_id


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Generate info YAML for use with single cell data QC pipeline

Usage:
  generate_info_yaml.py --md=<md> --outdir=<outdir> [options]

Arguments:
  -m --md=<md>              Metadata CSV file (required)
  -o --outdir=<outdir>      Output directory (required)

Options:
  -h --help                 Show this screen
"""


# ==============================
# FUNCTIONS
# ==============================
@logger.catch(reraise=True)
def _main(opt: dict) -> None:
    # Read input CSV and check fields are valid
    md = pd.read_csv(opt["--md"], header=0)
    assert set(md.columns).issuperset(
        {"donor", "pool", "sample", "hto", "cells_loaded"}
    ), "Invalid metadata CSV file."

    # Add unique sample ID
    md = md.assign(
        sample_id=lambda x: sample_id(x.donor.tolist(), x.pool.tolist()),
    )

    # Generate info YAML
    logger.info("Generating info YAML")
    generate_info_yaml(
        df=md[["sample_id", "hto", "sample", "cells_loaded"]],
        filename=os.path.join(opt["--outdir"], "info.yaml"),
    )
    logger.success("Output file: {}", os.path.abspath(os.path.join(opt["--outdir"], "info.yaml")))


def generate_info_yaml(df: pd.DataFrame, filename: str | None = None) -> dict:
    """
    Generate info YAML.

    Arguments:
        ``df``: DataFrame containing run metadata.\n
        ``filename``: Output file path or ``None``.

    Returns:
        Writes formatted YAML string to ``filename`` (if provided) and returns dictionary containing YAML data.
    """
    out = {
        sample_id: hto.loc[sample_id].to_dict("index")
        for sample_id, hto in df.set_index(["sample_id", "hto"]).groupby("sample_id")
    }

    if filename is not None:
        with open(file=filename, mode="w", encoding="UTF-8") as file:
            yaml.dump(data=out, stream=file)

    return out


# ==============================
# SCRIPT
# ==============================
if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
