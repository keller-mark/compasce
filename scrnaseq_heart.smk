include: "./common.smk"

ZARR_PATH = join(PROCESSED_DIR, "ht.adata.zarr")

# Rules
rule all:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.merged")

rule merge_metadata:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_pearson_residuals"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp"),
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance"),
    #join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_lemur")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.merged")
  resources:
    slurm_partition="short",
    runtime=60, # 1 hour
    mem_mb=4_000, # 4 GB
    cpus_per_task=2
  shell:
    """
    python scripts/merge_metadata.py \
        --zarr-path {ZARR_PATH}
    """

rule compute_lemur:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_lemur")
  resources:
    slurm_partition="medium",
    runtime=60*24*4, # 4 days
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_lemur" \
        --mem-limit "24GB" \
        --n-workers 4 \
        --threads-per-worker 2
    """

rule compute_diffexp:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffexp")
  resources:
    slurm_partition="medium",
    runtime=60*24*3, # 3 days
    mem_mb=160_000, # 160 GB
    cpus_per_task=4
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_diffexp"
    """

rule compute_diffabundance:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.compute_diffabundance")
  resources:
    slurm_partition="medium",
    runtime=60*24*2, # 2 days
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "compute_diffabundance"
    """

rule densmap:
  input:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.normalize_basic")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata.densmap")
  resources:
    slurm_partition="short",
    runtime=60*11, # 11 hours
    mem_mb=240_000, # 240 GB
    cpus_per_task=4
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
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
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
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    compasce \
        --zarr-path {ZARR_PATH} \
        --function-name "normalize_basic"
    """

rule convert_to_zarr:
  input:
    raw=join(RAW_DIR, "HT_raw.h5ad"),
    processed=join(RAW_DIR, "HT_processed.h5ad")
  output:
    join_zdone(ZARR_PATH, "uns", "comparison_metadata")
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    python scripts/run_comparisons_heart.py \
        --input-raw {input.raw} \
        --input-processed {input.processed} \
        --output {ZARR_PATH} \
        --stop-early \
        --no-subset
    """

rule download_proc_h5ad:
  output:
    join(RAW_DIR, "HT_processed.h5ad")
  params:
    file_url="https://hubmap-data-products.s3.amazonaws.com/4ea0c58e-13a9-41fd-b698-de9124c02873/HT_processed.h5ad"
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    curl -L -o {output} "{params.file_url}"
    """

rule download_raw_h5ad:
  output:
    join(RAW_DIR, "HT_raw.h5ad")
  params:
    file_url="https://hubmap-data-products.s3.amazonaws.com/4ea0c58e-13a9-41fd-b698-de9124c02873/HT_raw.h5ad"
  resources:
    slurm_partition="short",
    runtime=60*2, # 2 hours
    mem_mb=120_000, # 120 GB
    cpus_per_task=2
  shell:
    """
    curl -L -o {output} "{params.file_url}"
    """
