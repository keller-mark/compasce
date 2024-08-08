import pytest

import zarr
import numpy as np
import pandas as pd
import os
from os.path import join
from anndata import read_h5ad
import dask.array as da
import shutil

import comparisce as csc

from .fixtures import DATA_DIR, adata_fixture, client_fixture

@pytest.fixture
def writing_zarr_path():
    zarr_path = join(DATA_DIR, "test_zarr_writing.h5ad.zarr")
    shutil.rmtree(zarr_path, ignore_errors=True)
    return zarr_path

@pytest.fixture
def creation_zarr_path():
    zarr_path = join(DATA_DIR, "test_zarr_creation.h5ad.zarr")
    shutil.rmtree(zarr_path, ignore_errors=True)
    return zarr_path

def get_dirsize(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return total_size

@pytest.mark.limit_memory("30 MB")
def test_zarr_creation(adata_fixture, client_fixture, creation_zarr_path):
    client = client_fixture
    adata = adata_fixture
    zarr_path = creation_zarr_path

    csc.io.create_lazy_anndata(adata, zarr_path, client=client)
    assert os.path.exists(join(zarr_path, "layers", "counts", ".zdone"))
    assert get_dirsize(zarr_path) >= 60 * 1000000 # 60 MB

def test_zarr_writing(adata_fixture, client_fixture, writing_zarr_path):
    client = client_fixture
    adata = adata_fixture
    zarr_path = writing_zarr_path

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







    