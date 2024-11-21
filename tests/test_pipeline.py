import pytest

import zarr
import numpy as np
import pandas as pd
from os.path import join
from anndata import read_h5ad
import dask.array as da
import shutil

import compasce as csc

from .fixtures import DATA_DIR, adata_fixture, client_fixture

@pytest.fixture
def normalization_zarr_path():
    zarr_path = join(DATA_DIR, "test_normalization.h5ad.zarr")
    shutil.rmtree(zarr_path, ignore_errors=True)
    return zarr_path

@pytest.fixture
def diffexp_zarr_path():
    zarr_path = join(DATA_DIR, "test_diffexp.h5ad.zarr")
    shutil.rmtree(zarr_path, ignore_errors=True)
    return zarr_path

@pytest.mark.skip(reason="slow")
def test_normalization(adata_fixture, client_fixture, normalization_zarr_path):
    client = client_fixture
    adata = adata_fixture
    zarr_path = normalization_zarr_path
    ladata = csc.io.create_lazy_anndata(adata, zarr_path, client=client)

    # Normalize basic
    csc.normalize_basic(ladata)

    # Check the results
    z = zarr.open(zarr_path, mode="r")
    X_logcounts = z["/layers/logcounts"]

    assert X_logcounts.shape == (20000, 3000)
    assert X_logcounts.dtype in (np.float32, np.float64)
    assert X_logcounts.chunks == (20000, 5)
    assert np.sum(X_logcounts) == pytest.approx(27327000.0)

    # Normalize pearson residuals
    csc.normalize_pearson_residuals(ladata)

    # Check the results
    X_pearson_residuals = z["/layers/pearson_residuals"]

    assert X_pearson_residuals.shape == (20000, 3000)
    assert X_pearson_residuals.dtype in (np.float32, np.float64)
    assert X_pearson_residuals.chunks == (20000, 5)
    assert np.nansum(X_pearson_residuals) == pytest.approx(127542.29)

    # Run densMAP
    csc.densmap(ladata)

    X_densmap = z["/obsm/X_densmap"]
    assert X_densmap.shape == (20000, 2)
    assert X_densmap.dtype in (np.float32, np.float64)
    assert np.sum(X_densmap) == pytest.approx(112528.72)
    
def test_diffexp(adata_fixture, client_fixture, diffexp_zarr_path):
    client = client_fixture
    adata = adata_fixture
    zarr_path = diffexp_zarr_path

    cdata = csc.io.ComparativeData(zarr_path=zarr_path)
    ladata = cdata.create_lazy_anndata(adata, client=client)

    # Normalize basic
    csc.normalize_basic(ladata)

    # Compute diffexp
    csc.compute_diffexp(cdata, ladata)

    # Check the results
    z = zarr.open(zarr_path, mode="r")

    assert z['/compare_cell_type.val_monocyte.__rest__/ranked_genes.adata.zarr/var'].attrs == {
        "_index": "_index",
        "column-order": [
            "names",
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj"
        ],
        "encoding-type": "dataframe",
        "encoding-version": "0.2.0"
    }

    z_attrs_dict = dict(z.attrs)
    assert 'consolidated_uns' in z_attrs_dict
    assert 'compare_cell_type.val_b_cell.__rest__' in z_attrs_dict['consolidated_uns']
    assert 'ranked_genes.adata.zarr' in z_attrs_dict['consolidated_uns']['compare_cell_type.val_b_cell.__rest__']
    assert 'ranked_pathways.adata.zarr' in z_attrs_dict['consolidated_uns']['compare_cell_type.val_b_cell.__rest__']
    assert len(z_attrs_dict['consolidated_uns'].keys()) == 26

    assert z_attrs_dict['consolidated_uns']['compare_cell_type.val_b_cell.__rest__']['ranked_genes.adata.zarr'] == {
        "obsType": "cell",
        "featureType": "gene",
        "obsSetSelection": [["cell_type", "B cell"]],
    }

