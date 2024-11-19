# compasce

<!--[![PyPI](https://img.shields.io/pypi/v/compasce)](https://pypi.org/project/compasce)-->

Data processing functions to support visual comparisons of single-cell data.


:warning: Work in progress



## Installation

```sh
pip install compasce
```

## Usage

```python
import compasce as csc

adata = read_h5ad("my_adata.h5ad")
zarr_path = "my_adata.h5ad.zarr"
client = csc.create_dask_client()

csc.run_all(adata, zarr_path, client=client)
```

Or, run functions individually:

```python
import compasce as csc

adata = read_h5ad("my_adata.h5ad")
zarr_path = "my_adata.h5ad.zarr"
client = csc.create_dask_client()

ladata = csc.io.create_lazy_anndata(adata, zarr_path, client=client)

# Normalization
csc.normalize_basic(ladata)
csc.normalize_pearson_residuals(ladata)
```

## ComparativeData format

We define a ComparativeData object which is a container for AnnData, MuData, and SpatialData objects.
This format serves as a convention for how to organize pre-computed comparison results and store them on-disk.
It will not support all possible comparative use cases, but instead aims to support a set of common use cases that we have identified.

The comparative data format, again, is a container Zarr store for existing formats from the scverse ecosystem, themselves also stored via Zarr.

<!-- raw text for https://tree.nathanfriend.com
my_atlas.cdata.zarr
  - __all__                          # no comparison or filtering
    - cells.adata.zarr             # TODO: also support mudata
      - uns/comparasce           # special metadata, will also be zarr-consolidated
          obsType: "cell"
          featureType: "gene"
    - participants.adata.zarr    
      - uns/comparasce
        obsType: "participant"
        featureType: "clinical"    
  - compare_celltype.val_b_cell.__rest__
    - ranked_genes.adata.zarr
      - uns/comparasce
        featureType: "gene"
    - ranked_pathways.adata.zarr
      - uns/comparasce
        featureType: "pathway"
  - compare_celltype.val_cytotoxic_t_cell.__rest__
    - ranked_genes.adata.zarr
    - ranked_pathways.adata.zarr
  - compare_disease.val_healthy_reference.val_aki
     - adata.zarr                         # lemur results
  - filter_celltype.val_fibroblast.compare_disease.val_healthy_reference.val_AKI
    - ranked_genes.adata.zarr
    - ranked_pathways.adata.zarr
  - .zmetadata                              # zarr consolidated metadata
  - .zattrs                                 # uns-consolidated metadata
-->

```
my_atlas.cdata.zarr
├── __all__                          # no comparison or filtering
│   ├── cells.adata.zarr             # TODO: also support mudata
│   │   └── uns/comparasce           # special metadata, will also be zarr-consolidated
│   │       ├── obsType: "cell"
│   │       └── featureType: "gene"
│   └── participants.adata.zarr    
│       └── uns/comparasce
│           ├── obsType: "participant"
│           └── featureType: "clinical"    
├── compare_celltype.val_b_cell.__rest__
│   ├── ranked_genes.adata.zarr
│   │   └── uns/comparasce
│   │       └── featureType: "gene"
│   └── ranked_pathways.adata.zarr
│       └── uns/comparasce
│           └── featureType: "pathway"
├── compare_celltype.val_cytotoxic_t_cell.__rest__
│   ├── ranked_genes.adata.zarr
│   └── ranked_pathways.adata.zarr
├── compare_disease.val_healthy_reference.val_aki
│   └── adata.zarr                         # lemur results
├── filter_celltype.val_fibroblast.compare_disease.val_healthy_reference.val_AKI
│   ├── ranked_genes.adata.zarr
│   └── ranked_pathways.adata.zarr
├── .zmetadata                              # zarr consolidated metadata
└── .zattrs                                 # uns-consolidated metadata
```

Principles:
- The intermediate directory names should follow a formula, but are primarily meant to be human readable (as opposed to machine readable). Downstream apps/tools should not rely on these names.
- Machine-readable metadata should be stored in the `uns/compasce` dictionaries, which are then consolidated into `/.zattrs` and `/.zmetadata` in the root of the `.cdata.zarr` store.
- A downstream application should be able to read the `/.zmetadata` and `/.zattrs` data to understand all comparisons that were performed and where that data is stored within the rest of the zarr store.

### Note on alternative approaches

Other approaches seemed to have limitations, for example storage of a single differential expression test result under `adata.uns["rank_genes_groups"]` using numpy structured arrays or returning results in a standalone pandas dataframe.
These approaches work well when plots are generated using python and but we need more standard ways to organize such results in order to develop interactive tools around them.
Another alternative would be to port the comparative methods to webassembly or javascript but this space of methods moves rapidly and porting/compilation can be time consuming.
Some methods may also have long execution times or high computational resource requirements, prohibiting such an approach.




## Development

```sh
uv sync --extra dev
uv run pytest
```

Download example data:

```sh
bash ./scripts/download_data.sh
```

Create a subset of example data:

```sh
python ./scripts/subset_h5ad.py \
    --input ./data/lake_et_al.full.h5ad \
    --output ./data/lake_et_al.subset.h5ad
```

### Testing

```sh
pytest
```
