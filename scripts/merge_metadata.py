from compasce.io import MultiComparisonMetadata, LazyAnnData
from compasce import create_o2_dask_client
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr-path", type=str, required=True, help = "Path to zarr store.")
    args = parser.parse_args()

    zarr_path = args.zarr_path
    client = create_o2_dask_client(memory_limit='1GB')

    cm = MultiComparisonMetadata()
    cm.load_state(zarr_path, include_comparisons=True)
    cm.merge_states(zarr_path, suffixes=[
        'normalize_basic',
        'normalize_pearson_residuals',
        'densmap',
        'compute_diffexp',
        'compute_diffabundance',
        'compute_lemur',
    ])

    ladata = LazyAnnData(zarr_path, client=client)

    ladata.uns["comparison_metadata"] = cm.serialize()
    ladata.save(arr_path=["uns", "comparison_metadata.merged"])