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
    parser.add_argument("--overwrite", type=bool, default=False, required=False)
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

        # Cleanup of sample-level data
        # Reference samples are mapped to the empty string in the PrimaryAdjudicatedCategory column
        adata.obs["PrimaryAdjudicatedCategory"] = adata.obs["PrimaryAdjudicatedCategory"].apply(lambda v: "Reference" if v == "" else v)
        adata.obs["EnrollementCategoryAlt"] = adata.obs["EnrollementCategory"].apply(lambda v: "Reference" if v in ["LD", "HRT"] else v)
        return adata

    # For KPMP_PREMIERE....h5ad
    sample_id_col = "SampleID"
    sample_group_pairs = [
        # AKI vs. HRT
        ('EnrollementCategoryAlt', ('Reference', 'AKI')),
        # AKI vs. H-CKD
        ('EnrollementCategoryAlt', ('AKI', 'H-CKD')),
        # D-CKD vs. HRT
        ('EnrollementCategoryAlt', ('DKD', 'Reference')),
        ('PrimaryAdjudicatedCategory', ('Diabetic Kidney Disease', 'Reference')),
        # Acute tubular injury vs. HRT
        ('PrimaryAdjudicatedCategory', ('Acute Tubular Injury', 'Reference')),
        # Acute interstitial nephritis vs. HRT
        ('PrimaryAdjudicatedCategory', ('Acute Interstitial Nephritis', 'Reference')),
        # Diabetes CKD vs. Hypertension CKD
        ('EnrollementCategoryAlt', ('DKD', 'H-CKD')),
        ('PrimaryAdjudicatedCategory', ('Diabetic Kidney Disease', 'Hypertensive Kidney Disease')),
        # ATN vs. AIN
        ('PrimaryAdjudicatedCategory', ('Acute Interstitial Nephritis', 'Acute Tubular Injury')),
        # Broader categories
        ('diseasetype', ('Reference', 'AKI')),
        ('diseasetype', ('CKD', 'AKI')),
        ('diseasetype', ('Reference', 'CKD')),
    ]

    ladata = run_all(
        get_adata,
        zarr_path=args.output,
        overwrite=args.overwrite,
        client=create_o2_dask_client(memory_limit=args.mem_limit),
        sample_id_col=sample_id_col,
        sample_group_pairs=sample_group_pairs,
    )

    print("Done")
