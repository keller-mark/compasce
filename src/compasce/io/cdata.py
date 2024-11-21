import zarr
from os.path import join, dirname, basename
import re
import glob
import numpy as np

from .zarr_io import dispatched_write_zarr
from .lazy_anndata import LazyAnnData
from ..constants import COMPASCE_KEY


def slugify(text):
    text = text.lower()
    text = re.sub(r'[^\w\s-]', '', text)  # Remove non-word characters
    text = re.sub(r'[\s]+', '_', text)  # Replace spaces with underscores
    return text

def dir_name_to_str(dir_name):
    if isinstance(dir_name, str):
        assert dir_name in ["__all__"]
        return dir_name
    elif isinstance(dir_name, list):
        result = []
        for item in dir_name:
            if isinstance(item, tuple):
                assert len(item) == 2
                assert item[0] in ["filter", "compare", "val"]
                item_key = item[0]
                item_val = slugify(item[1])
                result.append(f"{item_key}_{item_val}")
            else:
                assert isinstance(item, str)
                assert item in ["__rest__", "__all__"]
                result.append(item)
        return ".".join(result)

def simplify_value(val):
    # AnnData stores list values in uns as numpy arrays,
    # which cause problems when zarr tries to serialize the dict
    # for the standard zarr zattrs.
    if isinstance(val, dict):
        return {k: simplify_value(v) for k, v in val.items()}
    elif isinstance(val, np.ndarray):
        return val.tolist()
    return val

class ComparativeData:
    def __init__(self, zarr_path, group_pairs=[]):
        self.zarr_path = zarr_path
        # Zarr open_consolidated will fail unless the root contains a group.
        zarr.open(zarr_path, mode="a")
        self.group_pairs = []
    
    def create_lazy_anndata(self, adata, dir_name="__all__", name=None, arr_path=None, mode="w", client=None):
        adata_path = join(self.zarr_path, dir_name_to_str(dir_name), f"{name}.adata.zarr" if name is not None else "adata.zarr")
        dispatched_write_zarr(adata, adata_path, arr_path=arr_path, mode=mode, client=client)

        return LazyAnnData(adata_path, client=client)

    def update(self):
        # Merge the dicts saved in /*/*.adata.zarr/uns/compasce, and store them in the root attrs
        uns_dict = {}
        adata_paths = glob.glob(join(self.zarr_path, "*", "*.adata.zarr"))
        for adata_path in adata_paths:
            adata_name = basename(adata_path)
            group_name = basename(dirname(adata_path))
            ladata = LazyAnnData(adata_path)
            if group_name not in uns_dict:
                uns_dict[group_name] = {}
            uns_dict[group_name][adata_name] = simplify_value(ladata.uns.get(COMPASCE_KEY))

        z = zarr.open(self.zarr_path, mode="r+")
        z.attrs["consolidated_uns"] = uns_dict

        # Conslidate zarr metadata
        zarr.consolidate_metadata(self.zarr_path)
        