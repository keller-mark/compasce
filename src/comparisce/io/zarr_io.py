import numpy as np
import zarr
import dask.array as da
from anndata import AnnData
from anndata._io.zarr import read_dataframe, _read_legacy_raw, _clean_uns
from anndata.experimental import read_dispatched, write_dispatched, read_elem


def dispatched_read_zarr(store, omit_X=False, omit_raw=True, omit_layers=None):
    # Function that reads an AnnData object from a Zarr store but omits certain keys.
    # Adapted from https://github.com/scverse/anndata/blob/1461fecd1712eefb1e5a5c0a75547b0e169a23d5/src/anndata/_io/zarr.py#L51
    if isinstance(store, zarr.Group):
        f = store
    else:
        f = zarr.open(store, mode="r")

    # Read with handling for backwards compat
    def callback(func, elem_name: str, elem, iospec):
        #print(f"Reading {elem_name}")
        if elem_name == "/X" and omit_X:
            return None
        if elem_name == "/raw" and omit_raw:
            return None
        if omit_layers is not None and elem_name in [f"/layers/{l}" for l in omit_layers]:
            return None

        if elem_name == "/layers":
            # We want to trick the anndata read_basic_zarr function into thinking that this Zarr group
            # only contains a subset of keys.
            # Reference: https://github.com/scverse/anndata/blob/1461fecd1712eefb1e5a5c0a75547b0e169a23d5/src/anndata/_io/specs/methods.py#L145
            elem = {
                k: v
                for k, v in elem.items()
                if omit_layers is None or k not in omit_layers
            }

        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            return AnnData(
                **{
                    k: read_dispatched(v, callback)
                    for k, v in elem.items()
                    if not k.startswith("raw.")
                }
            )
        elif elem_name.startswith("/raw."):
            return None
        elif elem_name in {"/obs", "/var"}:
            return read_dataframe(elem)
        elif elem_name == "/raw":
            # Backwards compat
            return _read_legacy_raw(f, func(elem), read_dataframe, func)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    # Backwards compat (should figure out which version)
    if "raw.X" in f:
        raw = AnnData(**_read_legacy_raw(f, adata.raw, read_dataframe, read_elem))
        raw.obs_names = adata.obs_names
        adata.raw = raw

    # Backwards compat for <0.7
    if isinstance(f["obs"], zarr.Array):
        _clean_uns(adata)

    return adata


def dispatched_write_zarr(adata, out_path, var_chunk_size=5, arr_path=None, mode="r+"):
    # Write to Zarr and set custom chunk shape for layers
    # Reference: https://anndata.readthedocs.io/en/latest/tutorials/notebooks/%7Bread%2Cwrite%7D_dispatched.html
    def write_chunked(func, store, k, elem, dataset_kwargs, iospec):
        """Write callback that chunks X and layers"""

        def set_chunks(d, chunks=None):
            """Helper function for setting dataset_kwargs. Makes a copy of d."""
            d = dict(d)
            if chunks is not None:
                d["chunks"] = chunks
            else:
                d.pop("chunks", None)
            return d

        if iospec.encoding_type == "array":
            if 'layers' in k or k.endswith('X'):
                dataset_kwargs = set_chunks(dataset_kwargs, (adata.shape[0], var_chunk_size))
            else:
                dataset_kwargs = set_chunks(dataset_kwargs, None)

        if isinstance(elem, da.Array):
            raise ValueError("Dask arrays should not be passed to write_chunked")
        elif elem is None:
            print("Skipping writing of None element")
        else:
            func(store, k, elem, dataset_kwargs=dataset_kwargs)

    z = zarr.open(out_path, mode=mode)

    # Monkey patch the clear method to prevent clearing the root group
    # Reference: https://github.com/scverse/anndata/blob/1461fecd1712eefb1e5a5c0a75547b0e169a23d5/src/anndata/_io/specs/registry.py#L299C1-L299C26
    old_clear = z.clear
    z.clear = (lambda: None) # Do not allow clearing the root group

    old_delitem = z.__class__.__delitem__
    def patched_delitem(self, item):
        print(f"Attepting to delete {item}")
        if item == "/layers":
            pass
        else:
            old_delitem(self, item)
    z.__class__.__delitem__ = patched_delitem

    write_dispatched(z, "/", adata, callback=write_chunked)
    # Restore (though not really necessary)
    z.clear = old_clear
    z.__class__.__delitem__ = old_delitem


    if arr_path is not None:
        if "write_metadata" not in z["uns"]:
            write_metadata = z["uns"].create_group("write_metadata")
        else:
            write_metadata = z["uns"]["write_metadata"]
        
        group_name = arr_path[0]
        subgroup_name = arr_path[1] if len(arr_path) > 1 else None
        if group_name not in write_metadata:
            group = write_metadata.create_group(group_name)
        else:
            group = write_metadata[group_name]
        
        if subgroup_name is not None:
            group.create_group(subgroup_name)    