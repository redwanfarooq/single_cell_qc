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

# Set module rules list
module_rules = ['droplet_qc', 'libraries_qc', 'ambient_removal', 'summary']

# Import rules
include: 'rules/droplet_qc.smk'
include: 'rules/libraries_qc.smk'
include: 'rules/ambient_removal.smk'
include: 'rules/summary.smk'

# Set targets list
targets = [x for rule in [droplet_qc, libraries_qc, ambient_removal, summary] for x in rule]
# --------------------------------------------------


# --------------------------------------------------
# RULES
rule all:
	input: [os.path.abspath(x) for x in targets]
# --------------------------------------------------