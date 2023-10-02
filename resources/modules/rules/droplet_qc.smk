##########################################################################################
# Snakemake rule for droplet processing QC
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

rmd_dir = config.get("rmd_dir", "resources/rmarkdown")

# Define rule
rule droplet_qc:
    output: os.path.join(config["output_dir"]["reports"], "droplet_qc/{sample}.html")
    log: os.path.abspath("logs/droplet_qc/{sample}.log")
    threads: 1
    params:
        path = rmd_dir if os.path.isabs(rmd_dir) else os.path.join(workflow.basedir, rmd_dir),
        metadata = lambda wildcards: get_hto_metadata(wildcards, info=info),
        gex_matrix = config["gex_matrix"],
        hto_matrix = config["hto_matrix"],
        output = os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}.csv"),
        emptydrops_lower = config.get("emptydrops_lower", None),
        emptydrops_niters = config.get("emptydrops_niters", None),
        demuxmix_model = config.get("demuxmix_model", None),
        demuxmix_pAcpt = config.get("demuxmix_pAcpt", None),
        multiplet_projection = config.get("multiplet_projection", None)
    conda: "singlecell-r"
    envmodules: "R-cbrg"
    message: "Running droplet processing QC for {wildcards.sample}"
    script: "{params.path}/droplet_qc.Rmd"

droplet_qc = expand("{path}/droplet_qc/{sample}.html", path=[config["output_dir"]["reports"]], sample=samples)
