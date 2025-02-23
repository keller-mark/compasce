from compasce import run_all, create_dask_client, create_o2_dask_client
from anndata import read_h5ad
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True, help = "Path to KPMP Premiere H5AD file.")
    parser.add_argument("--output", type=str, required=True, help = "Path to output zarr store directory")
    parser.add_argument("--subset", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--mem-limit", type=str, default='16GB', required=False)
    parser.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--stop-early", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()

    def get_adata():
        adata = read_h5ad(args.input)

        should_subset = args.subset
        if should_subset:
            print("SUBSETTING")
            # subset using random sample so that multiple sample groups are represented to enable comparison
            np.random.seed(1)
            obs_subset = np.random.choice(adata.obs.index.tolist(), size=20_000, replace=False).tolist()
            var_slice = slice(None)
            adata = adata[obs_subset, var_slice].copy()
            adata.layers["counts"] = adata.raw[obs_subset, var_slice].X.todense()
        else:
            print("NOT SUBSETTING")
            adata.layers["counts"] = adata.raw.X.todense()
        adata.raw = None

        # Cleanup of sample-level data
        # Reference samples are mapped to the empty string in the AdjudicatedCategory column
        adata.obs["AdjudicatedCategory"] = adata.obs["PrimaryAdjudicatedCategory"].apply(lambda v: "Reference" if v == "" else v)
        adata.obs["EnrollmentCategory"] = adata.obs["EnrollementCategory"].apply(lambda v: "Reference" if v in ["LD", "HRT"] else v)
        return adata

    # For KPMP_PREMIERE....h5ad
    donor_id_col = "donor_id"
    sample_id_col = "SampleID"
    sample_group_pairs = [
        # AKI vs. HRT
        ('EnrollmentCategory', ('Reference', 'AKI')),
        # AKI vs. H-CKD
        ('EnrollmentCategory', ('AKI', 'H-CKD')),
        # D-CKD vs. HRT
        ('EnrollmentCategory', ('DKD', 'Reference')),
        ('AdjudicatedCategory', ('Diabetic Kidney Disease', 'Reference')),
        # Acute tubular injury vs. HRT
        ('AdjudicatedCategory', ('Acute Tubular Injury', 'Reference')),
        # Acute interstitial nephritis vs. HRT
        ('AdjudicatedCategory', ('Acute Interstitial Nephritis', 'Reference')),
        # Diabetes CKD vs. Hypertension CKD
        ('EnrollmentCategory', ('DKD', 'H-CKD')),
        ('AdjudicatedCategory', ('Diabetic Kidney Disease', 'Hypertensive Kidney Disease')),
        # ATN vs. AIN
        ('AdjudicatedCategory', ('Acute Interstitial Nephritis', 'Acute Tubular Injury')),
        # Broader categories
        ('diseasetype', ('Reference', 'AKI')),
        ('diseasetype', ('CKD', 'AKI')),
        ('diseasetype', ('Reference', 'CKD')),
    ]
    cell_type_cols = [
        "cell_type",
        "subclass.l1",
        "subclass.l2",
    ]

    ladata = run_all(
        get_adata,
        zarr_path=args.output,
        overwrite=args.overwrite,
        client=create_o2_dask_client(memory_limit=args.mem_limit),
        donor_id_col=donor_id_col,
        sample_id_col=sample_id_col,
        sample_group_pairs=sample_group_pairs,
        cell_type_cols=cell_type_cols,
        stop_early=args.stop_early,
    )

    print("Done")
