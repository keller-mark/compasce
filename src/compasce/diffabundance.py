import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
#import pertpy as pt
from pertpy.tools._coda._sccoda import Sccoda


def compute_diffabundance(ladata, cm):
    print(f"Running diff abundance tests")

    # Check for a .zdone file
    if ladata.has_zdone(["uns", "compute_diffabundance"]):
        return ladata

    cell_type_col = cm.cell_type_col
    sample_group_pairs = cm.sample_group_pairs
    sample_id_col = cm.sample_id_col

    cell_types = ladata.obs[cell_type_col].unique().tolist()
    cell_types = [x for x in cell_types if pd.notna(x)]

    sccoda_model = Sccoda()

    for sample_group_pair in sample_group_pairs:
        sample_group_col, (sample_group_left, sample_group_right) = sample_group_pair
        print(f"Running scCODA for {sample_group_col}: {sample_group_left} vs. {sample_group_right}")

        cmp = cm.add_comparison([("compare", sample_group_col), ("val", sample_group_left), ("val", sample_group_right)])

        try:
            sccoda_data = sccoda_model.load(
                ladata,
                type="cell_level",
                generate_sample_level=True,
                cell_type_identifier=cell_type_col,
                sample_identifier=sample_id_col,
                covariate_obs=[sample_group_col],
            )

            # Select control and case data
            sccoda_data.mod["coda_case_vs_control"] = sccoda_data["coda"][
                sccoda_data["coda"].obs[sample_group_col].isin([sample_group_left, sample_group_right])
            ].copy()


            sccoda_data_2 = sccoda_model.prepare(
                sccoda_data,
                modality_key="coda_case_vs_control",
                formula=sample_group_col,
                reference_cell_type="automatic",
            )

            # Run MCMC
            sccoda_model.run_nuts(sccoda_data_2, modality_key="coda_case_vs_control")

            intercept_df = sccoda_model.get_intercept_df(sccoda_data_2, modality_key="coda_case_vs_control")
            effect_df = sccoda_model.get_effect_df(sccoda_data_2, modality_key="coda_case_vs_control")
            # Remove multi-index of effect_df
            effect_df = effect_df.reset_index(level='Covariate')

            joint_df = pd.merge(effect_df, intercept_df, on="Cell Type", suffixes=('_effect', '_intercept'), validate="one_to_one")
            
            assert joint_df.shape[0] == effect_df.shape[0] # Should have the same number of rows.

            sccoda_params = sccoda_data_2.mod["coda_case_vs_control"].uns["scCODA_params"]
            method_params = {
                "reference_cell_type": sccoda_params["reference_cell_type"],
                "automatic_reference_absence_threshold": sccoda_params["automatic_reference_absence_threshold"],
            }

            uns_key = cmp.append_df("uns", "sccoda_df", method_params, {
                "obsType": "cell",
                "sampleSetSelection": [[sample_group_col, sample_group_right]],
                "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
            })
            ladata.uns[uns_key] = joint_df

        except ValueError:
            print(f"Error while running scCODA for {sample_group_col}: {sample_group_left} vs. {sample_group_right}")
    
    ladata.write_zdone(["uns", "compute_diffabundance"])
    return ladata
