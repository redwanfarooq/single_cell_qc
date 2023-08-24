#!/bin/env python


"""
Runs single cell data QC pipeline.
"""


# ==============================
# MODULES
# ==============================
import os
import yaml
import docopt


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC = """
Run single cell data QC pipeline

Usage:
  run.py [options]

Options:
  -h --help                 Show this screen
  -u --update               Update module scripts using rule specifications in 'config/modules.yaml' (will not run pipeline)
"""


# ==============================
# GLOBAL VARIABLES
# ==============================
CONSOLE_LOG = "echo $(date '+[%Y-%m-%d%t%H:%M:%S]') {}"


# ==============================
# FUNCTIONS
# ==============================
def _main(opt: dict) -> None:
    # Get and execute shell command
    cmd = _get_cmd(update=opt["--update"])
    os.system(" && ".join(cmd))


def _cmd(*args, message: str):
    cmd = [" ".join(args)]
    cmd.insert(0, CONSOLE_LOG.format(message))
    cmd.append(CONSOLE_LOG.format("Done."))
    return cmd


def _get_cmd(update: bool = False) -> list[str]:
    if update:
        cmd = _cmd(
            f"{SCRIPTS_DIR}/generate_modules.py",
            "--modules=config/modules.yaml",
            "--template=resources/templates/module.template",
            "--outdir=resources/modules",
            message="Updating module scripts...",
        )
    else:
        cmd = _cmd(
            f"{SCRIPTS_DIR}/generate_wrapper.py",
            f"--module={MODULE}",
            "--template=resources/templates/wrapper.template",
            message="Generating wrapper script...",
        )
        cmd += _cmd(
            f"{SCRIPTS_DIR}/generate_info_yaml.py",
            f"--md={SAMPLES_CSV}",
            f"--outdir={METADATA_DIR}",
            message="Generating info YAML file...",
        )
        cmd += _cmd(
            "snakemake --profile=profile",
            message="Starting single cell data QC pipeline...",
        )
    return cmd


# ==============================
# SCRIPT
# ==============================
with open(file="config/config.yaml", mode="r", encoding="UTF-8") as file:
    config = yaml.load(stream=file, Loader=yaml.SafeLoader)
    SCRIPTS_DIR = config.get("scripts_dir", "resources/scripts")
    METADATA_DIR = config.get("metadata_dir", "metadata")
    try:
        SAMPLES_CSV = os.path.join(METADATA_DIR, config["samples"])
        MODULE = config["module"]
    except KeyError as err:
        raise KeyError(f"{err} not specified in '{file.name}'") from err


with open(file="config/modules.yaml", mode="r", encoding="UTF-8") as file:
    try:
        RULES = yaml.load(stream=file, Loader=yaml.SafeLoader)[MODULE]
    except KeyError as err:
        raise KeyError(f"Module {err} not specified in '{file.name}'") from err

if __name__ == "__main__":
    _main(opt=docopt.docopt(DOC))
