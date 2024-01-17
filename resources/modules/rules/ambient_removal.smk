##########################################################################################
# Snakemake rule for ambient removal
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/droplet_qc.smk
# Requires outputs from resources/rules/libraries_qc.smk
##########################################################################################

rmd_dir = config.get("rmd_dir", "resources/rmarkdown")

# Define rule
rule ambient_removal:
    input: os.path.join(config["output_dir"]["reports"], "libraries_qc/{sample}.html")
    output: os.path.join(config["output_dir"]["reports"], "ambient_removal/{sample}.html")
    log: os.path.abspath("logs/ambient_removal/{sample}.log")
    threads: 1
    params:
        path = rmd_dir if os.path.isabs(rmd_dir) else os.path.join(workflow.basedir, rmd_dir),
        droplet_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}"),
        libraries_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "libraries_qc/{sample}"),
        outdir = lambda wildcards: os.path.join(config["output_dir"]["data"], "ambient_removal/{sample}"),
        gex_matrix = lambda wildcards: config.get("gex_matrix", None),
        atac_matrix = lambda wildcards: config.get("atac_matrix", None),
        adt_matrix = lambda wildcards: config.get("adt_matrix", None),
        control = config.get("control", None),
        gex_ambient_removal = config.get("gex_ambient_removal", None),
        adt_ambient_removal = config.get("adt_ambient_removal", None)
    conda: "single_cell_qc"
    envmodules: "R-cbrg"
    message: "Running ambient removal for {wildcards.sample}"
    script: "{params.path}/ambient_removal.Rmd"

ambient_removal = expand("{path}/ambient_removal/{sample}.html", path=[config["output_dir"]["reports"]], sample=samples)