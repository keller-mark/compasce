from .normalization import normalize_basic, normalize_pearson_residuals
from .densmap import densmap
from .diffexp import compute_diffexp
from .io.cdata import ComparativeData
from .constants import COMPASCE_KEY


def run_all(get_adata, zarr_path, client=None, sample_id_col=None, sample_group_pairs=None):
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

    if "counts" not in adata.layers:
        raise ValueError("adata.layers['counts'] must exist")
    if adata.shape[0] < 1:
        raise ValueError("adata must have at least one row")
    if adata.shape[1] < 1:
        raise ValueError("adata must have at least one column")

    cdata = ComparativeData(
        zarr_path=zarr_path,
        sample_group_pairs=sample_group_pairs,
        sample_id_col=sample_id_col
    )
    
    adata.uns[COMPASCE_KEY] = {
        "obsType": "cell",
        "featureType": "gene",
        "featureValueType": "expression",
    }
    ladata = cdata.create_lazy_anndata(adata, client=client)
    cdata.create_sample_anndata(ladata)

    del adata

    # depends on: uns/write_metadata/layers/counts
    # creates: uns/write_metadata/layers/logcounts
    normalize_basic(ladata)

    del ladata.layers["counts"]

    # depends on: uns/write_metadata/layers/counts
    # creates: /layers/pearson_residuals
    normalize_pearson_residuals(ladata)

    del ladata.layers["logcounts"]

    densmap(ladata)

    # depends on: uns/write_metadata/layers/logcounts
    # creates: varm/DE_cell_type_vs_rest
    compute_diffexp(cdata, ladata)

    
    # for group_pair in group_pairs:
    #     compute_diffexp(ladata, group_pair)
    #     compute_lemur(cdata, group_pair)

    return True