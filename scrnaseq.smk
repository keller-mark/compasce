include: "./common.smk"

ZARR_PATH = join(PROCESSED_DIR, "kpmp_premiere.adata.zarr")

# Rules
rule all:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_pearson_residuals"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_lemur")

rule compute_lemur:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_lemur")
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_lemur"
    """

rule compute_diffabundance:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance")
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_diffabundance"
    """

rule compute_diffexp:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp")
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_diffexp"
    """

rule densmap:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap")
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "densmap"
    """

rule normalize_pearson_residuals:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_pearson_residuals")
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "normalize_pearson_residuals"
    """

rule normalize_basic:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "normalize_basic"
    """

rule convert_to_zarr:
  input:
    join(RAW_DIR, "KPMP_PREMIERE_SC_version1.5_ForExplorer_withRC.032624.h5ad")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  shell:
    """
    python scripts/run_comparisons.py \
        --input {input} \
        --output {ZARR_PATH} \
        --stop-early \
        --subset
    """

