# comparisce

[![PyPI](https://img.shields.io/pypi/v/comparisce)](https://pypi.org/project/comparisce)

Preprocessing functions that support comparison of single-cell data

## Installation

```sh
pip install comparisce
```

## Usage

```python
import comparisce as csc

adata = read_h5ad("my_adata.h5ad")
zarr_path = "my_adata.h5ad.zarr"
client = csc.create_dask_client()

csc.run_all(adata, zarr_path, client=client)
```

Or, run functions individually:

```python
import comparisce as csc

adata = read_h5ad("my_adata.h5ad")
zarr_path = "my_adata.h5ad.zarr"
client = csc.create_dask_client()

ladata = csc.io.create_lazy_anndata(adata, zarr_path, client=client)

# Normalization
csc.normalize_basic(ladata)
csc.normalize_pearson_residuals(ladata)
```

## ComparativeData format

A ComparativeData object is a container for AnnData, MuData, and SpatialData objects.
We define this format as a convention for how to store pre-computed comparison results.
We acknowledge that it does not support all possible comparative use cases, however we have designed it to support common use cases.

## Development

```sh
conda create -n comparisce-dev python=3.12
conda activate comparisce-dev
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
