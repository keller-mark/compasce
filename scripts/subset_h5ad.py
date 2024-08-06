import argparse
from anndata import read_h5ad
import numpy as np

NUM_CELLS = 20000
NUM_GENES = 3000

random_state = 1234

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    adata = read_h5ad(args.input)
    
    np.random.seed(random_state)
    random_cells = np.random.choice(adata.obs.index, NUM_CELLS, replace=False)
    random_genes = np.random.choice(adata.var.index, NUM_GENES, replace=False)

    adata.layers = {}
    adata.layers["counts"] = adata.raw.X.toarray()
    
    # Subset the AnnData object
    adata = adata[random_cells, random_genes].copy()

    adata.raw = None

    adata.write(args.output)

