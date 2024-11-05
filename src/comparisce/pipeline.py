import os
import platform

from .normalization import normalize_basic, normalize_pearson_residuals
from .densmap import densmap
from .diffexp import compute_diffexp
from .io.lazy_anndata import create_lazy_anndata
from .io.cdata import ComparativeData


def run_all(adata, zarr_path, client=None):
    """
    adata = read_h5ad("path/to/adata.h5ad")
    zarr_path = "path/to/adata.zarr"
    client = create_dask_client()

    run_all(adata, zarr_path, client=client)
    """

    cdata = ComparativeData(zarr_path=zarr_path)
    
    ladata = cdata.create_lazy_anndata(adata, client=client)

    # depends on: uns/write_metadata/layers/counts
    # creates: uns/write_metadata/layers/logcounts
    normalize_basic(ladata)

    # depends on: uns/write_metadata/layers/counts
    # creates: /layers/pearson_residuals
    normalize_pearson_residuals(ladata)

    densmap(ladata)

    # depends on: uns/write_metadata/layers/logcounts
    # creates: varm/DE_cell_type_vs_rest
    compute_diffexp(ladata)

    return True