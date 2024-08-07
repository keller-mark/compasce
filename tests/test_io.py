import pytest

import zarr
import numpy as np
import pandas as pd
from os.path import join
from anndata import read_h5ad
import dask.array as da
import shutil

import comparisce as csc

DATA_DIR = join("data")
ADATA_PATH = join(DATA_DIR, "lake_et_al.subset.h5ad")

def test_zarr_writing():
    adata = read_h5ad(ADATA_PATH)
    zarr_path = join(DATA_DIR, "test_zarr_writing.h5ad.zarr")
    shutil.rmtree(zarr_path, ignore_errors=True)

    client = csc.create_dask_client(memory_limit="2GB")

    ladata = csc.io.create_lazy_anndata(adata, zarr_path, client=client)
    ladata.save()

    # Check that things have been written properly.
    z = zarr.open(zarr_path, mode="r")

    assert str(z.tree()).strip().startswith("""
/
 ├── X
 │   ├── data (3593680,) float32
 │   ├── indices (3593680,) int32
 │   └── indptr (20001,) int32
 ├── layers
 │   └── counts (20000, 3000) float32
 ├── obs
 """.strip())

    # Check that things have been written properly.
    X_dask = ladata.get_da_from_zarr_layer("counts")

    assert X_dask.shape == (20000, 3000)
    assert X_dask.chunksize == (20000, 5)

    X_sum = da.sum(X_dask).compute()
    print(X_sum)
    assert X_sum == pytest.approx(13361987.0)

    X_double = X_dask * 2

    # Write dask array
    ladata.put_da_to_zarr_layer("counts_double", X_double)

    # Read using plain zarr
    X_sum_double = np.sum(z['/layers/counts_double'][()])

    assert X_sum_double == pytest.approx(26723974.0)







    