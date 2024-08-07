import scanpy as sc
import pandas as pd
import numpy as np

def densmap(ladata, random_state=1234):
    # Neighbors will use /X or /obsm/X_pca by default.

    # TODO: implement a mechanism via LazyAnnData to map a particular layer to X
    sc.pp.neighbors(ladata, random_state=random_state)
    sc.tl.densmap(ladata, random_state=random_state)

    ladata.save(arr_path=["obsm/X_densmap"])

    return ladata
