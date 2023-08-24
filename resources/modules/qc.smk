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
dirs = ['droplet_qc']
# --------------------------------------------------


# --------------------------------------------------
# RULES
rule all:
	input: expand("{path}/{dir}/{sample}.html", path=[config["output_dir"]["reports"]], dir=dirs, sample=samples)
# --------------------------------------------------