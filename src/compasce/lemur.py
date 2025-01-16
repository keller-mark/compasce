import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData

import pylemur

from .constants import COMPASCE_KEY


def compute_lemur(cdata, ladata):

    sample_group_pairs = cdata.sample_group_pairs

    def get_input_arr():
        return ladata.get_da_from_zarr_layer("logcounts")

    for sample_group_pair in sample_group_pairs:
        sample_group_col, (sample_group_left, sample_group_right) = sample_group_pair
        model = pylemur.tl.LEMUR(ladata, get_input_arr, design = f"~ {sample_group_col}", n_embedding=15, layer = "logcounts", copy=False)
        model.fit()
        model.align_with_harmony(max_iter=3) # TODO: remove this max_iter param

        ctrl_pred = model.predict(new_condition=model.cond(**{ sample_group_col: sample_group_left }))
        stim_pred = model.predict(new_condition=model.cond(**{ sample_group_col: sample_group_right }))
        lemur_diff_matrix = (stim_pred - ctrl_pred)

        # Recalculate the DensMAP on the embedding calculated by LEMUR
        ladata.obsm["lemur_embedding"] = model.embedding
        ladata.save(arr_path=["obsm", "lemur_embedding"])

        # Map a particular layer to X, to trick the scanpy functions that only work with X.
        ladata.set_alias(on_disk_path=["obsm", "lemur_embedding"], in_mem_path=["X"])

        sc.pp.neighbors(ladata, n_neighbors=30, use_rep="X", key_added="lemur_embedding_neighbors")
        sc.tl.umap(ladata, method="densmap", key_added="lemur_densmap", neighbors_key="lemur_embedding_neighbors")

        # Store results in a new AnnData object.
        # TODO: copy over obs/var index columns?
        lemur_adata = AnnData(X=lemur_diff_matrix, obs=None, var=None, obsm={"X_densmap": ladata.obsm["lemur_densmap"]})
        lemur_adata.uns[COMPASCE_KEY] = {
            "obsType": "cell",
            "featureType": "gene",
            "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
        }
        cdata.create_lazy_anndata(lemur_adata, dir_name=[("compare", sample_group_col), ("val", sample_group_left), ("val", sample_group_right)], name="lemur")

    return ladata
