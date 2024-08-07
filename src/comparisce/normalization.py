import scanpy as sc
import numpy as np
import dask.array as da

from .dask import create_dask_wrapper


def normalize_basic(ladata, input_key="counts", output_key="logcounts"):
    ladata.layers[output_key] = ladata.layers[input_key].copy()

    # Scanpy gets confused by non-AnnData objects even when it is AnnData-like
    adata = ladata.adata
    sc.pp.normalize_total(adata, target_sum = 1e6, layer=output_key, inplace=True)
    sc.pp.log1p(adata, layer=output_key, copy=False)

    ladata.save(arr_path=["layers", output_key])

    return ladata


def _normalize_pearson_residuals(X, theta = 100, clip = None):
    # check theta
    if theta <= 0:
        # TODO: would "underdispersion" with negative theta make sense?
        # then only theta=0 were undefined..
        raise ValueError('Pearson residuals require theta > 0')

    # prepare clipping
    if clip is None:
        n = X.shape[0]
        clip = da.sqrt(n)
    if clip < 0:
        raise ValueError("Pearson residuals require `clip>=0` or `clip=None`.")

    sums_genes = da.sum(X, axis=0, keepdims=True)
    sums_cells = da.sum(X, axis=1, keepdims=True)
    sum_total = da.sum(sums_genes)

    mu = da.divide(da.matmul(sums_cells, sums_genes), sum_total)
    diff = da.subtract(X, mu)
    residuals = da.divide(diff, da.sqrt(da.add(mu, da.divide(da.power(mu, 2), theta))))
    return residuals.clip(min=-clip, max=clip)


def normalize_pearson_residuals(adata, input_key="counts", output_key="pearson_residuals"):
    def get_input_arr():
        return adata.get_da_from_zarr_layer(input_key)
    
    def put_output_arr(output_arr):
        adata.put_da_to_zarr_layer(output_key, output_arr)

    normalize_pearson_residuals_dask = create_dask_wrapper(_normalize_pearson_residuals)
    normalize_pearson_residuals_dask(get_input_arr, put_output_arr)

    return adata
