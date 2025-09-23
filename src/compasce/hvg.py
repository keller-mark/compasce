import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
import pylemur


def annotate_hvgs(ladata, cm):

    if ladata.has_zdone(["uns", "annotate_hvgs"]):
        return ladata

    hvg_df = sc.pp.highly_variable_genes(
        ladata,
        layer="logcounts",
        n_top_genes=1000,
        flavor='seurat',
        # TODO: should a batch-specific key be used?
        batch_key=cm.sample_id_col,
        subset=False,
        inplace=False
    )
    hvg_cols = [
        'highly_variable',
        'means',
        'dispersions',
        'dispersions_norm',
        'variances',
        'variances_norm',
        'highly_variable_rank',
        'highly_variable_nbatches',
        'highly_variable_intersection',
    ]
    for colname in hvg_cols:
        if colname in hvg_df.columns:
            ladata.var[colname] = hvg_df[colname]
    
    ladata.obsm["logcounts_hvg"] = ladata.layers["logcounts"].get_orthogonal_selection((
        # All cells
        slice(None),
        # Subset to the highly variable genes using boolean mask
        hvg_df['highly_variable'].values
    ))

    ladata.save(arr_path=["obsm", "X_hvg"])
    ladata.write_zdone(["uns", "annotate_hvgs"])

    return ladata
