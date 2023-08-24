##########################################################################################
# Snakemake rule for droplet processing QC
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule droplet_qc:
    output: os.path.join(config["output_dir"]["reports"], "droplet_qc/{sample}.html")
    log: os.path.abspath("logs/droplet_qc/{sample}.log")
    threads: 1
    params:
        path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
        metadata = lambda wildcards: get_hto_metadata(wildcards, info=info),
        gex_dir = config["gex_dir"],
        hto_dir = config["hto_dir"],
        output_dir = os.path.join(config["output_dir"]["data"], "droplet_qc"),
        emptydrops_lower = config.get("emptydrops_lower", 100),
        emptydrops_niters = config.get("emptydrops_niters", 10000),
        demuxmix_model = config.get("demuxmix_model", "auto"),
        demuxmix_pAcpt = config.get("demuxmix_pAcpt", 0.9)
    conda: "sc"
    envmodules: "R-cbrg"
    message: "Running droplet processing QC for {wildcards.sample}"
    script: "{params.path}/rmd/droplet_qc.Rmd"
