from compasce import run_all, create_dask_client, create_o2_dask_client
from anndata import read_h5ad
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True, help = "Path to KPMP Premiere H5AD file.")
    parser.add_argument("--output", type=str, required=True, help = "Path to output zarr store directory")
    parser.add_argument("--subset", type=bool, default=False, required=False)
    parser.add_argument("--mem-limit", type=str, default='16GB', required=False)
    args = parser.parse_args()

    def get_adata():
        adata = read_h5ad(args.input)

        should_subset = args.subset
        if should_subset:
            # subset using random sample so that multiple sample groups are represented to enable comparison
            np.random.seed(1)
            obs_subset = np.random.choice(adata.obs.index.tolist(), size=20_000, replace=False).tolist()
            var_slice = slice(None)
            adata = adata[obs_subset, var_slice].copy()
            adata.layers["counts"] = adata.raw[obs_subset, var_slice].X.todense()
        else:
            adata.layers["counts"] = adata.raw.X.todense()
        adata.raw = None
        return adata

    # For KPMP_PREMIERE....h5ad
    sample_id_col = "SampleID"
    sample_group_pairs = [
        ('diseasetype', ('Reference', 'AKI')),
        ('diseasetype', ('CKD', 'AKI')),
        ('diseasetype', ('Reference', 'CKD')),
        # TODO: expand? use adjudicated?
    ]

    ladata = run_all(
        get_adata,
        zarr_path=args.output,
        client=create_o2_dask_client(memory_limit=args.mem_limit),
        sample_id_col=sample_id_col,
        sample_group_pairs=sample_group_pairs,
    )

    print("Done")
