import pytest

import zarr
import numpy as np
import pandas as pd
from os.path import join
from anndata import read_h5ad
import dask.array as da
import shutil

import compasce as csc

DATA_DIR = join("data")
ADATA_PATH = join(DATA_DIR, "lake_et_al.subset.h5ad")

@pytest.fixture(scope="module")
def client_fixture():
    return csc.create_dask_client(memory_limit="2GB")

@pytest.fixture(scope="module")
def adata_fixture():
    return read_h5ad(ADATA_PATH)

