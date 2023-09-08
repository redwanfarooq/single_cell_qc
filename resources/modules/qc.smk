##########################################################################################
# Snakemake module
# Author: Redwan Farooq
# Module name: qc
##########################################################################################


# --------------------------------------------------
# SETUP
# Load modules
import os
import yaml
from resources.scripts.rule import *

# Load and parse sample/library info from YAML file
with open(file=os.path.join(config.get("metadata_dir", "metadata"), "info.yaml"), mode="r", encoding="UTF-8") as file:
    info = yaml.load(stream=file, Loader=yaml.SafeLoader)
for key, value in parse_info(info).items():
    globals()[key] = value

# Import rules
include: 'rules/droplet_qc.smk'
include: 'rules/libraries_qc.smk'
include: 'rules/batch_correction.smk'

# Set targets
targets = [x for rule in [droplet_qc, libraries_qc, batch_correction] for x in rule]
# --------------------------------------------------


# --------------------------------------------------
# RULES
rule all:
	input: [os.path.abspath(x) for x in targets]
# --------------------------------------------------