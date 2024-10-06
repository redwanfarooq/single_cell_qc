##########################################################################################
# Snakemake rule for droplet processing QC
# Author: Redwan Farooq
# Requires functions from resources/scripts/rule.py
##########################################################################################

scripts_dir = config.get("scripts_dir", "resources/scripts")

# Define rule
rule droplet_qc:
    input: os.path.join(config["output_dir"]["reports"], "cellbender/{sample}.html") if "cellbender" in module_rules else os.path.join(config["output_dir"]["data"], f"merge/{{sample}}/{'filtered' if config.get('pre_filtered', None) else 'raw'}_feature_bc_matrix.h5")
    output: os.path.join(config["output_dir"]["reports"], "droplet_qc/{sample}.html")
    log: os.path.abspath("logs/droplet_qc/{sample}.log")
    threads: 1
    params:
        script_path = scripts_dir if os.path.isabs(scripts_dir) else os.path.join(workflow.basedir, scripts_dir),
        metadata = os.path.abspath(os.path.join(config.get("metadata_dir", "metadata"), config["samples"])),
        features_matrix = lambda wildcards: get_features_matrix(wildcards, data_dir=config["output_dir"]["data"], cellbender="cellbender" in module_rules, filtered=config.get("pre_filtered", None)),
        hto_matrix = lambda wildcards: config.get("hto_matrix", None),
        hto_prefix = config.get("hto_prefix", None),
        outdir = lambda wildcards: os.path.join(config["output_dir"]["data"], "droplet_qc/{sample}"), # DO NOT CHANGE - downstream rules will search for output files in this directory
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
    # envmodules: "R-cbrg"
    message: "Running droplet processing QC for {wildcards.sample}"
    script: "{params.script_path}/droplet_qc.Rmd"

droplet_qc = expand(os.path.join(config["output_dir"]["reports"], "droplet_qc/{sample}.html"), sample=samples)
