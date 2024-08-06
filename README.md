# comparisce

[![PyPI](https://img.shields.io/pypi/v/comparisce)](https://pypi.org/project/comparisce)

Preprocessing functions that support comparison of single-cell data

## Installation

```sh
pip install comparisce
```

## Usage

Usage follows the scverse API conventions.
Parameter names follow the R implementation of `miQC`.

```python
import comparisce as csc

# ...

csc.calculate_miqc(adata)
csc.filter_cells(adata)
```


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
