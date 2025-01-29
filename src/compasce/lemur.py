import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
import pylemur


def compute_lemur(ladata, cm):

    sample_group_pairs = cm.sample_group_pairs

    def get_input_arr():
        return ladata.get_da_from_zarr_layer("logcounts")

    for sample_group_pair in sample_group_pairs:
        sample_group_col, (sample_group_left, sample_group_right) = sample_group_pair
        cmp = cm.add_comparison([("compare", sample_group_col), ("val", sample_group_left), ("val", sample_group_right)])

        model = pylemur.tl.LEMUR(ladata, get_input_arr, design = f"~ {sample_group_col}", n_embedding=15, layer = "logcounts", copy=False)
        model.fit()
        model.align_with_harmony()

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

        method_params = {
            "design": f"~ {sample_group_col}",
            "n_embedding": 15,
            "layer": "logcounts"
        }
        # Store the results
        layer_key = cmp.append_df("layers", "lemur_diff", method_params, {
            "obsType": "cell",
            "featureType": "gene",
            "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
        })
        ladata.layers[layer_key] = lemur_diff_matrix

        obsm_key = cmp.append_df("obsm", "lemur_embedding", method_params, {
            "obsType": "cell",
            "featureType": "gene",
            "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
        })
        ladata.obsm[obsm_key] = model.embedding

        obsm_key = cmp.append_df("obsm", "lemur_densmap", method_params, {
            "obsType": "cell",
            "featureType": "gene",
            "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
        })
        ladata.obsm[obsm_key] = ladata.obsm["lemur_densmap"]

    return ladata
