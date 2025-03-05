from compasce.io import MultiComparisonMetadata, write_zdone
import zarr
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr-path", type=str, required=True, help = "Path to zarr store.")
    args = parser.parse_args()

    zarr_path = args.zarr_path

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

    z = zarr.open(zarr_path, mode="a")
    z["/uns/comparison_metadata"] = cm.serialize()
    z["/uns/comparison_metadata"].attrs["encoding-type"] = "string"
    z["/uns/comparison_metadata"].attrs["encoding-version"] = "0.2.0"
    write_zdone(zarr_path, arr_path=["uns", "comparison_metadata.merged"])
