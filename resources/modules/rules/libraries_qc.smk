##########################################################################################
# Snakemake rule for libraries QC
# Author: Redwan Farooq
# Requires outputs from resources/rules/droplet_qc.smk
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule libraries_qc:
    input: os.path.join(config["output_dir"]["reports"], "droplet_qc/{sample}.html")
    output: os.path.join(config["output_dir"]["reports"], "libraries_qc/{sample}.html")
    log: os.path.abspath("logs/libraries_qc/{sample}.log")
    threads: 1
    params:
        script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
        droplet_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}"),
        outdir = lambda wildcards: os.path.join(config["output_dir"]["data"], "libraries_qc/{sample}"),
        gex_matrix = lambda wildcards: config.get("gex_matrix", None),
        atac_matrix = lambda wildcards: config.get("atac_matrix", None),
        adt_matrix = lambda wildcards: config.get("adt_matrix", None),
        gex_filters = config.get("gex_filters", None),
        atac_filters = config.get("atac_filters", None),
        adt_filters = config.get("adt_filters", None)
    conda: "single_cell_qc"
    # envmodules: "R-cbrg"
    message: "Running libraries QC for {wildcards.sample}"
    script: "{params.script_path}/libraries_qc.Rmd"

libraries_qc = expand("{path}/libraries_qc/{sample}.html", path=[config["output_dir"]["reports"]], sample=samples)
