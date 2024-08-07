import scanpy as sc
import pandas as pd
import numpy as np

def densmap(ladata, random_state=1234):
    # Neighbors will use /X or /obsm/X_pca by default.

    ladata.refresh(needs_X=True, needs_layers=[])
    adata = ladata.adata

    # TODO: implement a mechanism via LazyAnnData to map a particular layer to X
    sc.pp.neighbors(adata, random_state=random_state)
    sc.tl.densmap(adata, random_state=random_state)

    ladata.save(arr_path=["obsm/X_densmap"])

    return ladata
