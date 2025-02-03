import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
from requests.exceptions import ConnectionError

from pertpy.tools._enrichment import Enrichment
import blitzgsea as blitz

def compute_diffexp(ladata, cm):
    print(f"Running diffexp tests for cell types vs rest")

    # Check for a .zdone file
    if ladata.has_zdone(["uns", "compute_diffexp"]):
        return ladata

    cell_type_col = cm.cell_type_col

    key_added = "rank_genes_groups"
    sc.tl.rank_genes_groups(ladata, groupby=cell_type_col, method="wilcoxon", layer="logcounts", key_added=key_added)

    cell_types = ladata.obs[cell_type_col].unique().tolist()
    cell_types = [x for x in cell_types if pd.notna(x)]

    # Reference: https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/enrichment.html#Using-custom-gene-sets
    # Ensure that the gene nomenclature in your target sets is compatible with your .var_names.
    targets = blitz.enrichr.get_library("GO_Molecular_Function_2021")
    pt_enricher = Enrichment()
    enrichment_dict = pt_enricher.hypergeometric(ladata, targets=targets)
    #pt_enricher.score(ladata, targets=targets, layer="logcounts", key_added="pertpy_enrichment")
    # Returns a dictionary with clusters as keys and data frames of test results sorted on q-value as the items.
    #enrichment = pt_enricher.gsea(ladata, targets=targets, key_added="pertpy_enrichment_gsea")


    for cell_type in cell_types:
        print(f"Getting diffexp test results for {cell_type} vs rest")
        cmp = cm.add_comparison([("compare", cell_type_col), ("val", cell_type), "__rest__"])

        df = sc.get.rank_genes_groups_df(ladata, group=cell_type, key=key_added)
        df = df.sort_values(by="pvals_adj", ascending=False)

        uns_key = cmp.append_df("uns", "rank_genes_groups", {
            "rank_genes_groups": ladata.uns[key_added]["params"],
            "rank_genes_groups_df": {
                "group": cell_type,   
            },
        }, {
            "obsType": "cell",
            "featureType": "gene",
            "obsSetSelection": [[cell_type_col, cell_type]],
        })
        ladata.uns[uns_key] = df

        # Compute gene enrichment.
        enrichment_df = enrichment_dict[cell_type]
        uns_key = cmp.append_df("uns", "pertpy_hypergeometric", {
            "rank_genes_groups": ladata.uns[key_added]["params"],
            "pertpy_hypergeometric": {
                "group": cell_type,
                "pvals_adj_thresh": .05,
                "direction": "both",
                "corr_method": "benjamini-hochberg",
            },
        }, {
            "obsType": "cell",
            "featureType": "pathway",
            "obsSetSelection": [[cell_type_col, cell_type]],
        })
        ladata.uns[uns_key] = enrichment_df
        
        try:
            # Run sc.queries.enrich on the results.
            enrichment_df = sc.queries.enrich(ladata, group=cell_type, log2fc_min=2, pval_cutoff=.01)
            enrichment_df = enrichment_df.drop(columns=["query", "parents"])
            
            uns_key = cmp.append_df("uns", "enrich", {
                "rank_genes_groups": ladata.uns[key_added]["params"],
                "enrich": {
                    "group": cell_type,
                    "log2fc_min": 2,
                    "pval_cutoff": .01
                },
            }, {
                "obsType": "cell",
                "featureType": "pathway",
                "obsSetSelection": [[cell_type_col, cell_type]],
            })
            ladata.uns[uns_key] = enrichment_df
        except (AssertionError, ConnectionError):
            print(f"Gene enrichment query failed for cell type {cell_type}")
    
    del ladata.uns[key_added]
    

    # Within cell type (case vs. control)
    sample_group_pairs = cm.sample_group_pairs

    for cell_type in cell_types:
        print(f"Running diffexp test for {cell_type} and sample group pairs")
        for sample_group_pair in sample_group_pairs:

            sample_group_col, (sample_group_left, sample_group_right) = sample_group_pair
            cmp = cm.add_comparison([("filter", cell_type_col), ("val", cell_type), ("compare", sample_group_col), ("val", sample_group_left), ("val", sample_group_right)])
            try:
                ladata.obs["cell_type_sample_group"] = ladata.obs[cell_type_col].astype(str) + "_" + ladata.obs[sample_group_col].astype(str)
                sc.tl.rank_genes_groups(ladata, groupby="cell_type_sample_group", groups=[f"{cell_type}_{sample_group_right}"], reference=f"{cell_type}_{sample_group_left}", method="wilcoxon", layer="logcounts", key_added=key_added)

                df = sc.get.rank_genes_groups_df(ladata, group=f"{cell_type}_{sample_group_right}", key=key_added)
                df = df.sort_values(by="pvals_adj", ascending=False)

                uns_key = cmp.append_df("uns", "rank_genes_groups", {
                    "rank_genes_groups": ladata.uns[key_added]["params"],
                    "rank_genes_groups_df": {
                        "group": f"{cell_type}_{sample_group_right}",   
                    },
                }, {
                    "obsType": "cell",
                    "featureType": "gene",
                    "obsSetFilter": [[cell_type_col, cell_type]],
                    "sampleSetSelection": [[sample_group_col, sample_group_right]],
                    "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
                })
                ladata.uns[uns_key] = df

                # Compute gene enrichment.
                pt_enricher = Enrichment()
                enrichment_dict = pt_enricher.hypergeometric(ladata, targets=targets)

                enrichment_df = enrichment_dict[f"{cell_type}_{sample_group_right}"]
                uns_key = cmp.append_df("uns", "pertpy_hypergeometric", {
                    "rank_genes_groups": ladata.uns[key_added]["params"],
                    "pertpy_hypergeometric": {
                        "group": f"{cell_type}_{sample_group_right}",
                        "pvals_adj_thresh": .05,
                        "direction": "both",
                        "corr_method": "benjamini-hochberg",
                    },
                }, {
                    "obsType": "cell",
                    "featureType": "pathway",
                    "obsSetFilter": [[cell_type_col, cell_type]],
                    "sampleSetSelection": [[sample_group_col, sample_group_right]],
                    "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
                })
                ladata.uns[uns_key] = enrichment_df

                try:
                    # Run sc.queries.enrich on the results.
                    enrichment_df = sc.queries.enrich(ladata, group=f"{cell_type}_{sample_group_right}", log2fc_min=2, pval_cutoff=.01, key=key_added)
                    enrichment_df = enrichment_df.drop(columns=["query", "parents"])
                    
                    uns_key = cmp.append_df("uns", "enrich", {
                        "rank_genes_groups": ladata.uns[key_added]["params"],
                        "enrich": {
                            "group": f"{cell_type}_{sample_group_right}",
                            "log2fc_min": 2,
                            "pval_cutoff": .01
                        },
                    }, {
                        "obsType": "cell",
                        "featureType": "pathway",
                        "obsSetFilter": [[cell_type_col, cell_type]],
                        "sampleSetSelection": [[sample_group_col, sample_group_right]],
                        "sampleSetFilter": [[sample_group_col, sample_group_left], [sample_group_col, sample_group_right]],
                    })
                    ladata.uns[uns_key] = enrichment_df
                except (AssertionError, ConnectionError):
                    print(f"Gene enrichment query failed for cell type {cell_type} and sample group pair {sample_group_pair}")

                del ladata.uns[key_added]
            except (IndexError, ValueError) as e:
                print(f"Error: likely due to insufficient data for comparison for {cell_type} and sample group pair {sample_group_pair}")

    # TODO: within cell type (inside spatial region vs. outside)

    # Write a .zdone file
    ladata.write_zdone(["uns", "compute_diffexp"])

    return ladata
