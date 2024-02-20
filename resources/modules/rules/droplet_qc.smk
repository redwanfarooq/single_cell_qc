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
        gex_matrix = lambda wildcards: config["gex_matrix"],
        atac_matrix = lambda wildcards: config.get("atac_matrix", None),
        adt_matrix = lambda wildcards: config.get("adt_matrix", None),
        hto_matrix = lambda wildcards: config.get("hto_matrix", None),
        outdir = lambda wildcards: os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}"),
        pre_filtered = config.get("pre_filtered", None),
        cell_calling_algorithm = config.get("cell_calling_algorithm", None),
        emptydrops_umi_min = config.get("emptydrops_umi_min", None),
        emptydrops_umi_min_frac_median = config.get("emptydrops_umi_min_frac_median", None),
        emptydrops_cand_max_n = config.get("emptydrops_cand_max_n", None),
        emptydrops_ind_min = config.get("emptydrops_ind_min", None),
        emptydrops_ind_max = config.get("emptydrops_ind_max", None),
        emptydrops_niters = config.get("emptydrops_niters", None),
        multimodal_modalities = config.get("multimodal_modalities", None),
        multimodal_ordmag_quantile = config.get("multimodal_ordmag_quantile", None),
        multimodal_ordmag_ratio = config.get("multimodal_ordmag_ratio", None),
        demuxmix_model = config.get("demuxmix_model", None),
        demuxmix_pAcpt = config.get("demuxmix_pAcpt", None),
        multiplet_calling = config.get("multiplet_calling", None),
        scdblfinder_clusters = config.get("scdblfinder_clusters", None),
        scdblfinder_dbr_sd = config.get("scdblfinder_dbr_sd", None),
        multiplet_projection = config.get("multiplet_projection", None)
    conda: "single_cell_qc"
    envmodules: "R-cbrg"
    message: "Running droplet processing QC for {wildcards.sample}"
    script: "{params.path}/droplet_qc.Rmd"

droplet_qc = expand("{path}/droplet_qc/{sample}.html", path=[config["output_dir"]["reports"]], sample=samples)
