from .normalization import normalize_basic, normalize_pearson_residuals
from .densmap import densmap
from .diffexp import compute_diffexp
from .diffabundance import compute_diffabundance
from .lemur import compute_lemur
from .io.lazy_anndata import create_lazy_anndata, create_sample_df
from .io.comparison_metadata import MultiComparisonMetadata


def run_all(get_adata, zarr_path, client=None, sample_id_col=None, sample_group_pairs=None, cell_type_col="cell_type"):
    """
    def get_adata():
        return read_h5ad("path/to/adata.h5ad")
    zarr_path = "path/to/adata.zarr"
    client = create_dask_client()

    run_all(get_adata, zarr_path, client=client)
    """

    # We ask the caller to provide a function which returns `adata`, 
    # so that we can free the memory within the function.
    adata = get_adata()
    cm = MultiComparisonMetadata(
        sample_group_pairs=sample_group_pairs,
        sample_id_col=sample_id_col,
        cell_type_col=cell_type_col,
    )

    all_cmp = cm.add_comparison("__all__")
    # all_cmp.append_df("uns", "cells", None, { "obsType": "cell" }) # assumed

    if "counts" not in adata.layers:
        raise ValueError("adata.layers['counts'] must exist")
    if adata.shape[0] < 1:
        raise ValueError("adata must have at least one row")
    if adata.shape[1] < 1:
        raise ValueError("adata must have at least one column")

    ladata = create_lazy_anndata(adata, zarr_path, client=client)

    del adata

    sample_df = create_sample_df(ladata, cm)
    uns_key = all_cmp.append_df("uns", "samples", None, {
        "obsType": "sample"
    })
    ladata.uns[uns_key] = sample_df
    

    # depends on: uns/write_metadata/layers/counts
    # creates: uns/write_metadata/layers/logcounts
    normalize_basic(ladata)

    del ladata.layers["counts"]

    # depends on: uns/write_metadata/layers/counts
    # creates: /layers/pearson_residuals
    #normalize_pearson_residuals(ladata)

    del ladata.layers["logcounts"]

    densmap(ladata)

    compute_diffexp(ladata, cm)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()

    compute_diffabundance(ladata, cm)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()

    print("Starting normalize_pearson_residuals")
    normalize_pearson_residuals(ladata)
    ladata.save()

    print("Starting compute_lemur")
    compute_lemur(ladata, cm)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save()

    return ladata