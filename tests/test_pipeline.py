import pytest

import zarr
import numpy as np
import pandas as pd
from os.path import join
from anndata import read_h5ad
import dask.array as da

import comparisce as csc

DATA_DIR = join("data")
ADATA_PATH = join(DATA_DIR, "lake_et_al.subset.h5ad")

def test_normalization():
    adata = read_h5ad(ADATA_PATH)
    zarr_path = join(DATA_DIR, "test_normalization.h5ad.zarr")

    client = csc.create_dask_client(memory_limit="2GB")

    ladata = csc.io.create_lazy_anndata(adata, zarr_path, client=client)
    ladata.save()

    # Normalize basic
    csc.normalize_basic(ladata)

    # Check the results
    z = zarr.open(zarr_path, mode="r")
    X_logcounts = z["/layers/logcounts"]

    assert X_logcounts.shape == (20000, 3000)
    assert X_logcounts.dtype == np.float32
    assert X_logcounts.chunks == (20000, 5)
    assert np.sum(X_logcounts) == pytest.approx(27327000.0)

    # Normalize pearson residuals
    csc.normalize_pearson_residuals(ladata)

    # Check the results
    X_pearson_residuals = z["/layers/pearson_residuals"]

    assert X_pearson_residuals.shape == (20000, 3000)
    assert X_pearson_residuals.dtype == np.float32
    assert X_pearson_residuals.chunks == (20000, 5)
    assert np.nansum(X_pearson_residuals) == pytest.approx(127542.086)


    # csc.densmap(ladata)
    



