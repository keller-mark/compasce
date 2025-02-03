import scanpy as sc
import pandas as pd
import numpy as np

def densmap(ladata, random_state=1234):

    if ladata.has_zdone(["obsm", "X_densmap"]):
        return ladata

    # Map a particular layer to X, to trick the scanpy functions that only work with X.
    ladata.set_alias(on_disk_path=["layers", "logcounts"], in_mem_path=["X"])
    
    # Neighbors will choose between usage of /X or /obsm/X_pca unless use_rep is used.
    sc.pp.neighbors(ladata, random_state=random_state, use_rep="X", n_neighbors=30)
    sc.tl.umap(ladata, method="densmap", random_state=random_state)

    ladata.save(arr_path=["obsm", "X_densmap"])
    ladata.clear_aliases()

    return ladata
