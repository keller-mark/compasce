import zarr
from os.path import join
import re

from .zarr_io import dispatched_write_zarr
from .lazy_anndata import LazyAnnData


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

class ComparativeData:
    def __init__(self, zarr_path, group_pairs=[]):
        self.zarr_path = zarr_path
        self.group_pairs = []
    
    def create_lazy_anndata(self, adata, dir_name="__all__", name=None, arr_path=None, mode="w", client=None):
        adata_path = join(self.zarr_path, dir_name_to_str(dir_name), f"{name}.adata.zarr" if name is not None else "adata.zarr")
        dispatched_write_zarr(adata, adata_path, arr_path=arr_path, mode=mode, client=client)

        return LazyAnnData(adata_path, client=client)

    def update(self):
        try:
            zarr.consolidate_metadata(self.zarr_path)
        except zarr.errors.PathNotFoundError:
            pass
        # TODO: merge consolidated metadata with the dicts saved in /*/*.adata.zarr/uns/*
