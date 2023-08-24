# Description
Snakemake pipeline for single cell data QC
- Modularised workflow can be modified and/or extended for different experiment designs
- Add as a submodule in a bioinformatics project GitHub repository
```
git submodule add https://github.com/redwanfarooq/single_cell_qc single_cell_qc
```
- Update submodule to the latest version
```
git submodule update --remote
```

# Required software
1. Global environment
    - [Snakemake >=v7.31](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    - [docopt >v=0.6](https://github.com/docopt/docopt)
    - [pandas >v=2.0](https://pandas.pydata.org/docs/getting_started/install.html)
2. Specific modules
    - [R >=v4.2](https://cran.r-project.org)

# Setup
1. Install software for global environment (requires Anaconda or Miniconda - see [installation instructions](https://conda.io/projects/conda/en/stable/user-guide/install/index.html))
    - Download [environment YAML](/resources/envs/snakemake.yaml)
    - Create new conda environment from YAML
    ```
    conda env create -f snakemake.yaml
    ```
2. Install software for specific module(s)
    - Manually install required software from source and check that executables are available in **PATH** (using `which`) *and/or*
    - Create new conda environments with required software from YAML (as above - download [environment YAMLs](/resources/envs)) *and/or*
    - Check that required software is available to load as environment modules (using `module avail`)
3. Set up pipeline configuration file **config/config.yaml** (see comments in file for detailed instructions)
4. Set up profile configuration file **profile/config.yaml** (see comments in file for detailed instructions)

# Run
1. Activate global environment
```
conda activate snakemake
```
2. Execute **run.py** in root directory

# Input
Pipeline requires the following input files/folders:

## General

**REQUIRED:**

1. Preprocessed sequencing data as feature-barcode count matrices (either sparse matrix or CSV format)
2. Sample hashing summary table in CSV format with the following required fields (with headers):
- **donor**: donor ID
- **hash**: hash ID
- **sample**: sample name
- **hto**: hashtag antibody ID
- **cells_loaded**: number of cells loaded

# Output
Output directories will be created in specified locations with subfolders containing the output of each QC step specified in the module:
- **data**: processed data
- **reports:** HTML reports summarising QC metrics

# Modules

## Available modules
- qc

## Adding new module
1. Add entry to module rule specifications file **config/modules.yaml** with module name and list of rule names
2. Add additional rule definition files in **modules/rules** folder (if needed)
- Rule definition file **must** specify output file **Rule/Sample.html** in the **reports** output directory
- Rule definition file **must** have the same file name as the rule with the file extension **.smk**
3. Execute **run.py** in root directory with `--update` flag (needs to be repeated if there are any further changes to the module rule specification in **config/modules.yaml**)