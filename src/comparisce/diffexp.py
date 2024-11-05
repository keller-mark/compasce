import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData


# TODO: way to specify pairs of groups at sample-level for comparison and their corresponding column
def compute_diffexp(cdata, ladata, cell_type_col="cell_type"):
    key_added = "rank_genes_groups"
    sc.tl.rank_genes_groups(ladata, groupby=cell_type_col, method="wilcoxon", layer="logcounts", key_added=key_added)

    # TODO: run sc.queries.enrich on the results?

    cell_types = ladata.obs[cell_type_col].unique().tolist()
    for cell_type in cell_types:
        df = sc.get.rank_genes_groups_df(ladata, group=cell_type, key=key_added)
        df = df.sort_values(by="pvals_adj", ascending=False)

        # Store results in a new AnnData object's var dataframe.
        de_adata = AnnData(X=None, obs=None, var=df)
        de_adata.uns[key_added] = {
            "params": ladata.uns[key_added]["params"],
            "column": cell_type_col,
            "value": cell_type,
        }
        cdata.create_lazy_anndata(de_adata, dir_name=[("compare", cell_type_col), ("val", cell_type), "__rest__"], name="ranked_genes")

        # Compute gene enrichment.
        try:
            enrichment_df = sc.queries.enrich(ladata, group=cell_type, log2fc_min=2, pval_cutoff=.01)
            enrichment_df = enrichment_df.drop(columns=["query", "parents"])
            
            # Store results in a new AnnData object's var dataframe.
            enrichment_adata = AnnData(X=None, obs=None, var=enrichment_df)
            enrichment_adata.uns[key_added + "_enrich"] = {
                "params": ladata.uns[key_added]["params"],
                "column": cell_type_col,
                "value": cell_type,
            }
            cdata.create_lazy_anndata(enrichment_adata, dir_name=[("compare", cell_type_col), ("val", cell_type), "__rest__"], name="ranked_pathways")
        except AssertionError:
            print(f"Gene enrichment query failed for cell type {cell_type}")
    
    del ladata.uns[key_added]
    

    # TODO: within cell type (case vs control)

    cdata.update()

    return ladata

def compute_diffexp_all():
    pass