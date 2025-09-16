from compasce import run_all, create_dask_client, create_o2_dask_client
from anndata import read_h5ad
import numpy as np
import pandas as pd
import argparse
import h5py


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=str, required=True, help = "Path to KPMP H5AD file from Globus, august 2025.")
    parser.add_argument("--input-csv", type=str, required=True, help = "Path to KPMP clinical data CSV file.")
    parser.add_argument("--output", type=str, required=True, help = "Path to output zarr store directory")
    parser.add_argument("--subset", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--mem-limit", type=str, default='16GB', required=False)
    parser.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--stop-early", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()



    def get_adata():
        adata = read_h5ad(args.input_h5ad)

        should_subset = args.subset
        if should_subset:
            print("SUBSETTING")
            # subset using random sample so that multiple sample groups are represented to enable comparison
            np.random.seed(1)
            obs_subset = np.random.choice(adata.obs.index.tolist(), size=20_000, replace=False).tolist()
            var_slice = slice(None)
            adata = adata[obs_subset, var_slice].copy()
        else:
            print("NOT SUBSETTING")
        
        adata.layers["counts"] = adata.layers["counts"].todense()

        # Join adata.obs with clinical data from CSV
        clinical_data = pd.read_csv(args.input_csv)
        adata.obs = adata.obs.merge(clinical_data, left_on="patient", right_on="Participant ID", how="left")

        # Cleanup of sample-level data
        def clean_adjudicated_category(row):
            if row["Primary Adjudicated Category"] != "":
                return row["Primary Adjudicated Category"]
            else:
                # The row was empty, so perhaps this sample has not yet been adjudicated.
                # However, we also need to check that this was not a "Healthy Reference" sample,
                # as these never go through the adjudication process.
                if row["Enrollment Category"] in ["Healthy Reference"]:
                    return "Healthy Reference"
                return ""
        adata.obs["AdjudicatedCategory"] = adata.obs.apply(clean_adjudicated_category, axis='columns')
        adata.obs["EnrollmentCategory"] = adata.obs["Enrollment Category"]

        # TODO: process other clinical columns? Sex, age group, etc.

        adata.obs = adata.obs.rename(columns={"subclass.l1": "subclass_l1", "subclass.l2": "subclass_l2", "subclass.l3": "subclass_l3"})
        return adata

    donor_id_col = "patient"
    sample_id_col = "specimen"
    sample_group_pairs = [
        # AKI vs. HRT
        ('EnrollmentCategory', ('Healthy Reference', 'AKI')),
        # AKI vs. H-CKD. (H-CKD not in enrollment category values anymore. Should I use "Hypertension History" Yes/No column?)
        ('EnrollmentCategory', ('AKI', 'CKD')),
        # D-CKD vs. HRT. (D-CKD not in enrollment category values anymore. Should I use "Diabetes History" Yes/No column?)
        ('EnrollmentCategory', ('CKD', 'Healthy Reference')),
        # Diabetes CKD vs. Hypertension CKD. (DKD nor H-CKD not in enrollment category values anymore. Should I use Yes/No columns?)
        #('EnrollmentCategory', ('DKD', 'H-CKD')),
        # D-CKD vs. HRT
        ('AdjudicatedCategory', ('Diabetic Kidney Disease', 'Healthy Reference')),
        # Acute tubular injury vs. HRT
        ('AdjudicatedCategory', ('Acute Tubular Injury', 'Healthy Reference')),
        # Acute interstitial nephritis vs. HRT
        ('AdjudicatedCategory', ('Acute Interstitial Nephritis', 'Healthy Reference')),
        # Diabetes CKD vs. Hypertension CKD
        ('AdjudicatedCategory', ('Diabetic Kidney Disease', 'Hypertensive Kidney Disease')),
        # ATN vs. AIN
        ('AdjudicatedCategory', ('Acute Interstitial Nephritis', 'Acute Tubular Injury')),

        # TODO: use Diabetes History and Hypertension History columns here.
    ]
    cell_type_cols = [
        "subclass_l1",
        "subclass_l2",
        "subclass_l3",
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
