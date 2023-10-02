##########################################################################################
# Snakemake rule for batch correction
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
# Requires outputs from resources/rules/droplet_qc.smk
# Requires outputs from resources/rules/libraries_qc.smk
##########################################################################################

rmd_dir = config.get("rmd_dir", "resources/rmarkdown")

# Define rule
rule batch_correction:
    input: expand("{path}/{dir}/{sample}.html", path=[config["output_dir"]["reports"]], dir=["droplet_qc", "libraries_qc"], sample=samples)
    output: os.path.join(config["output_dir"]["reports"], "batch_correction/report.html")
    log: os.path.abspath("logs/batch_correction/log.log")
    threads: 1
    params:
        path = rmd_dir if os.path.isabs(rmd_dir) else os.path.join(workflow.basedir, rmd_dir),
        sample_id = samples,
        droplet_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}.csv"),
        libraries_qc = lambda wildcards: os.path.join(config["output_dir"]["data"], "libraries_qc/{sample}.csv"),
        output_dir = os.path.join(config["output_dir"]["data"], "batch_correction"),
        gex_matrix = lambda wildcards: config.get("gex_matrix", None),
        adt_matrix = lambda wildcards: config.get("adt_matrix", None),
        adt_filters = config.get("adt_filters", None),
        control = config.get("control", None),
        correction_method = config.get("correction_method", None),
        correction_projection = config.get("correction_projection", None)
    conda: "singlecell-r"
    envmodules: "R-cbrg"
    message: "Running batch correction"
    script: "{params.path}/batch_correction.Rmd"

batch_correction = [os.path.join(config["output_dir"]["reports"], "batch_correction/report.html")]
