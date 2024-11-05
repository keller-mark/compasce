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
This format serves as a convention for how to organize pre-computed comparison results.

It will not support all possible comparative use cases, but instead aims to support a set of common use cases that we have identified.



## Development

```sh
conda create -n compasce-dev python=3.12
conda activate compasce-dev
pip install -e ".[dev]"
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
