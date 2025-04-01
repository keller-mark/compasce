from compasce import run_all, create_dask_client, create_o2_dask_client
from anndata import read_h5ad
import numpy as np
import pandas as pd
import argparse
import h5py


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-raw", type=str, required=True, help = "Path to HuBMAP Heart data product H5AD file.")
    parser.add_argument("--input-processed", type=str, required=True, help = "Path to HuBMAP Heart data product H5AD file.")
    parser.add_argument("--output", type=str, required=True, help = "Path to output zarr store directory")
    parser.add_argument("--subset", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--mem-limit", type=str, default='16GB', required=False)
    parser.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--stop-early", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()

    # For HT_raw and HT_processed files
    donor_id_col = ""
    sample_id_col = "hubmap_id"
    sample_group_pairs = [
        ('sex', ('Male', 'Female')),
        ('race', ('White', 'Black or African American')),
        ('AgeGroup', ('50 and Above', 'Below 50')),
        ('BmiGroup', ('Underweight', 'Healthy Weight')),
        ('BmiGroup', ('Underweight', 'Overweight')),
        ('BmiGroup', ('Underweight', 'Obesity')),
        ('BmiGroup', ('Healthy Weight', 'Overweight')),
        ('BmiGroup', ('Healthy Weight', 'Obesity')),
        ('BmiGroup', ('Overweight', 'Obesity')),
        ('BmiGroup2', ('Underweight or Healthy Weight', 'Overweight or Obesity')),
    ]
    cell_type_cols = [
        "leiden",
        "predicted_label",
        "azimuth_label",
    ]

    def check_pair_has_enough_data(adata, sample_group_pair):
        [colname, value_pair] = sample_group_pair
        [left_value, right_value] = value_pair
        left_num_cells = adata.obs[adata.obs[colname] == left_value].shape[0]
        right_num_cells = adata.obs[adata.obs[colname] == right_value].shape[0]

        if left_num_cells == 0:
            print(f"Zero cells for value {left_value} in column {colname}")
        if right_num_cells == 0:
            print(f"Zero cells for value {right_value} in column {colname}")

    def verify_pair_has_enough_data(adata, sample_group_pair):
        [colname, value_pair] = sample_group_pair
        [left_value, right_value] = value_pair
        left_num_cells = adata.obs[adata.obs[colname] == left_value].shape[0]
        right_num_cells = adata.obs[adata.obs[colname] == right_value].shape[0]

        assert left_num_cells > 0, f"Zero cells for value {left_value} in column {colname}"
        assert right_num_cells > 0, f"Zero cells for value {right_value} in column {colname}"
    

    def get_adata():
        adata_raw = read_h5ad(args.input_raw)
        adata_processed = read_h5ad(args.input_processed)

        # filter the raw to the processed cells; then append the counts matrix to the processed anndata object
        cells_in_processed = adata_processed.obs.index.tolist()
        cell_mask = adata_raw.obs.index.isin(cells_in_processed)

        genes_in_processed = adata_processed.var.index.tolist()
        gene_mask = adata_raw.var.index.isin(genes_in_processed)

        filtered_adata_raw = adata_raw[cell_mask, gene_mask].copy()
        adata_processed.layers["counts"] = filtered_adata_raw.X.todense()
        adata = adata_processed

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

        # Fix categorical columns
        for colname in ["azimuth_id", "azimuth_label", "barcode", "cause_of_death", "cl_match_type", "dataset", "hubmap_id", "leiden", "predicted_CLID", "predicted_label", "race", "sex"]:
            adata.obs[colname] = adata.obs[colname].astype(str)
        
        # Cleanup of sample-level data
        def group_ages(v):
            if pd.notna(v):
                if v >= 50:
                    return "50 and Above"
                elif v < 50:
                    return "Below 50"
            return v

        # Reference: https://www.cdc.gov/bmi/adult-calculator/bmi-categories.html
        def group_bmis(v):
            if pd.notna(v):
                if v < 18.5:
                    return "Underweight"
                elif v >= 18.5 and v < 25:
                    return "Healthy Weight"
                elif v >= 25 and v < 30:
                    return "Overweight"
                elif v >= 30:
                    return "Obesity"
            return v

        adata.obs["AgeGroup"] = adata.obs["age"].apply(group_ages)
        adata.obs["BmiGroup"] = adata.obs["bmi"].apply(group_bmis)

        def group_bmis2(v):
            if pd.notna(v):
                if v < 25:
                    return "Underweight or Healthy Weight"
                elif v >= 25:
                    return "Overweight or Obesity"
            return v
        adata.obs["BmiGroup2"] = adata.obs["bmi"].apply(group_bmis2)

        for sample_group_pair in sample_group_pairs:
            check_pair_has_enough_data(adata, sample_group_pair)
        for sample_group_pair in sample_group_pairs:
            verify_pair_has_enough_data(adata, sample_group_pair)

        return adata

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
