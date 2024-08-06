import os
import platform

from .normalization import normalize_basic, normalize_pearson_residuals
from .io.lazy_anndata import create_lazy_anndata


def run_all(adata, zarr_path, client=None):
    """
    adata = read_h5ad("path/to/adata.h5ad")
    zarr_path = "path/to/adata.zarr"
    client = create_dask_client()

    run_all(adata, zarr_path, client=client)
    """

    ladata = create_lazy_anndata(adata, zarr_path, client=client)

    # depends on: uns/write_metadata/layers/counts
    # creates: uns/write_metadata/layers/logcounts
    normalize_basic(ladata)

    # depends on: uns/write_metadata/layers/counts
    # creates: /layers/pearson_residuals
    normalize_pearson_residuals(ladata)


    return True