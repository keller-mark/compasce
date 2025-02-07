# compasce

<!--[![PyPI](https://img.shields.io/pypi/v/compasce)](https://pypi.org/project/compasce)-->

**Compa**re **s**ingle **c**ell **e**xperiments: Example of a data processing pipeline which supports visualizations of case-vs-control-style comparisons of single-cell atlas data.


<!--
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

-->



<!--
## ComparativeData format

We define a ComparativeData object which is a container for AnnData, MuData, and SpatialData objects.
This format serves as a convention for how to organize pre-computed comparison results and store them on-disk.
It will not support all possible comparative use cases, but instead aims to support a set of common use cases that we have identified.

The ComparativeData object, on-disk, is a container Zarr store for existing formats from the scverse ecosystem, which themselves can also be stored via Zarr, leveraging the hierarchy of groups/arrays concepts.

**Note**: this format is subject to change as we gain experience using it downstream for visualization purposes or through other feedback.
-->
<!-- raw text for https://tree.nathanfriend.com
my_atlas.cdata.zarr
  - __all__                          # no comparison or filtering
    - cells.adata.zarr             # TODO: also support mudata
      - uns/compasce             # special metadata, will be uns-consolidated
          obsType: "cell"
          featureType: "gene"
    - participants.adata.zarr    
      - uns/compasce
        obsType: "participant"
        featureType: "clinical"    
  - compare_celltype.val_b_cell.__rest__
    - ranked_genes.adata.zarr
      - uns/compasce
        featureType: "gene"
    - ranked_pathways.adata.zarr
      - uns/compasce
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

<!--
```
my_atlas.cdata.zarr
├── __all__                          # no comparison or filtering
│   ├── cells.adata.zarr             # TODO: also support mudata
│   │   └── uns/compasce             # special metadata, will be uns-consolidated
│   │       ├── obsType: "cell"
│   │       └── featureType: "gene"
│   └── participants.adata.zarr    
│       └── uns/compasce
│           ├── obsType: "participant"
│           └── featureType: "clinical"    
├── compare_celltype.val_b_cell.__rest__
│   ├── ranked_genes.adata.zarr
│   │   └── uns/compasce
│   │       └── featureType: "gene"
│   └── ranked_pathways.adata.zarr
│       └── uns/compasce
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
- The intermediate directory names should follow a formula, but are primarily meant to be human readable (as opposed to machine readable). Downstream apps/tools should not rely on these names (but may rely on the names of the leaf `*.zarr` sub-sub-directories).
- Machine-readable metadata should be stored in the `uns/compasce` dictionaries, which are then consolidated into `/.zattrs` in the root of the `.cdata.zarr` store.
- A downstream application should be able to read the `/.zmetadata` and `/.zattrs` data to understand all comparisons that were performed and where that data is stored within the rest of the zarr store.

### Note on alternative approaches

Other approaches have limitations, for example, while a single differential expression test result can be stored in `adata.uns["rank_genes_groups"]` using numpy structured arrays, there is not a standard way to extend this to multiple test results or align the dataframe with other `var` metadata.
These approaches work well when plots are generated using python and but we need more standard ways to organize such results in order to develop interactive tools around them.
Another alternative would be to port the comparative methods to webassembly or javascript but this space of methods moves very rapidly and porting/compilation is often not trivial.
For example, methods may have long execution times or high computational resource requirements, complicating a porting approach.

-->

## Development

```sh
uv venv
source .venv/bin/activate
uv sync --extra dev
```

or

```sh
conda create -n compasce-env python=3.10
conda activate compasce-env
pip install -e .
```

<!--
Download example data:

```sh
curl -o ./data/lake_et_al.full.h5ad "https://datasets.cellxgene.cziscience.com/4ad02dd3-9773-49ac-a8fd-7e0581703b7d.h5ad"
```

Create a subset of example data:

```sh
python ./scripts/subset_h5ad.py \
    --input ./data/lake_et_al.full.h5ad \
    --output ./data/lake_et_al.subset.h5ad
```
-->

## Run

### HMS O2 cluster SLURM setup

```sh
srun -p interactive --pty -t 11:00:00 -n 4 --mem 164G bash
# source ~/.bashrc_mark
# ssh-add
conda activate compasce-env
```

### General setup

Here, we set a `DATA_DIR` variable which is useful when running the subsequent commands in a cluster environment.
Change this value depending on where you would like to store data files.

```sh
export DATA_DIR=/n/data1/hms/dbmi/gehlenborg/lab/scmd-analysis
```

Store H5AD file at `$DATA_DIR/raw`.
If input file is in `.h5Seurat` format, convert to `.h5ad` (for example, using [SeuratDisk](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)).

<!--
```sh
cd $DATA_DIR/raw
curl -L -o KPMP_PREMIERE_SC_version1.5_ForExplorer_withRC.032624.h5ad "https://storage.googleapis.com/vitessce-demo-data/kpmp-jan-2025/KPMP_PREMIERE_SC_version1.5_ForExplorer_withRC.032624.h5ad"
cd -
```
-->

<!--

```sh
uv run python scripts/run_comparisons.py \
  --input data/raw/KPMP_PREMIERE_SC_version1.5_ForExplorer_withRC.032624.h5ad \
  --output data/processed/kpmp_premiere_small.adata.zarr \
  --subset True \
  --overwrite True \
  --mem-limit 2GB
```
-->

```sh
python scripts/run_comparisons.py \
  --input $DATA_DIR/raw/KPMP_PREMIERE_SC_version1.5_ForExplorer_withRC.032624.h5ad \
  --output $DATA_DIR/processed/kpmp_premiere.adata.zarr \
  --subset False
```

See the `scripts/run_comparisons.py` script for additional command-line options.

If using `uv` to manage the python environment, prepend `uv run ` to the above command.

This script took approximately 48 hours to complete with 160 GB of RAM. It is not yet optimized to run independent pipeline steps in parallel.

<!--
### Upload the results

```sh
gcloud auth login

cd $DATA_DIR/processed

gsutil -m cp -r ./kpmp_premiere.adata.zarr gs://vitessce-demo-data/scmd-analysis-october-2023
```
-->

### Testing

```sh
uv run pytest
```
