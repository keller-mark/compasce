import zarr
import dask.array as da
from anndata import AnnData
from dask.distributed import progress

from anndata._core.aligned_mapping import LayersView

from .zarr_io import dispatched_read_zarr, dispatched_write_zarr


class MappingWrapper:

    def __init__(self, mapping):
        self.mapping = mapping

    def __getitem__(self, key):
        # TODO: return zarr array that 
        return self.mapping[key]
    
    def __setitem__(self, key, value):
        self.mapping[key][:, :] = value

class LazyAnnData(AnnData):

    zarr_path = None
    adata = None
    var_chunk_size = None
    client = None
    done_init = False
    z = None

    def __init__(self, zarr_path, var_chunk_size=5, client=None):
        self.zarr_path = zarr_path
        self.adata = dispatched_read_zarr(zarr_path)
        self.z = zarr.open(zarr_path, mode="r+")
        self.var_chunk_size = var_chunk_size
        self.client = client

        super().__init__()

        self.done_init = True
    
    def __getattribute__(self, key):
        orig_getattr = object.__getattribute__
        z = orig_getattr(self, "z")
        done_init = orig_getattr(self, "done_init")
        if not done_init:
            return orig_getattr(self, key)
        # __getattr__ only gets called for attributes that don't actually exist
        print(f"Getting {key}")
        if key in { "X" }:
            return z["/X"]
        elif key in { "obsm", "varm", "layers" }:
            subkeys = z[key].keys()
            print(f"Getting {key} subkeys: {subkeys}")
            return MappingWrapper({ subkey: z[f"/{key}/{subkey}"] for subkey in subkeys })
        elif key in {"obs", "var", "uns"}:
            return getattr(self.adata, key)
        
        return orig_getattr(self, key)
    
    def __setattr__(self, key, value):
        if not self.done_init:
            return super().__setattr__(key, value)
        print(f"Setting {key}")
        if key in {"obs", "var", "uns"}:
            setattr(self.adata, key, value)
        elif key in { "obsm", "varm", "layers", "X" }:
            raise ValueError(f"Cannot set {key} via LazyAnnData. Set on the zarr array directly.")
        
        super().__setattr__(key, value)
    

    def save(self, arr_path=None, mode="r+"):
        dispatched_write_zarr(self.adata, self.zarr_path, var_chunk_size=self.var_chunk_size, arr_path=arr_path, mode=mode, client=self.client)

    def copy_layer(self, src_key, dest_key):
        zarr.copy(self.z[f"/layers/{src_key}"], self.z["/layers"], name=dest_key)

    def get_da_from_zarr_layer(self, layer_key):
        X_path = f"/layers/{layer_key}"

        if self.client is None:
            raise ValueError("A Dask client must be provided to use get_da_from_zarr_layer")
            # TODO: instead of throwing, just run da.from_array() on the numpy array normally

        z = zarr.open(self.zarr_path, path=X_path, mode='r')

        # Use Zarr so that X does not get loaded into memory
        # which would cause Dask to include it in the task graph,
        # which causes the task graph to be too large.
        X_dask = da.from_zarr(url=self.zarr_path, component=X_path, chunks=z.chunks, inline_array=False)
        return X_dask

    def put_da_to_zarr_layer(self, layer_key, X_dask):

        X_path = f"/layers/{layer_key}"

        if self.client is None:
            raise ValueError("A Dask client must be provided to use put_da_to_zarr_layer")
            # TODO: instead of throwing, just run .compute() and write the numpy array normally

        # Store using Zarr in a distributed way so that we don't need to transfer data back from each worker,
        # which is too large for a single worker to handle.
        residuals_delayed = X_dask.to_zarr(url=self.zarr_path, component=X_path, overwrite=True, compute=False)
        residuals_future = self.client.persist(residuals_delayed)

        progress(residuals_future)

def create_lazy_anndata(adata, zarr_path, client=None):
    assert "counts" in adata.layers, "The AnnData object must have a 'counts' layer."

    dispatched_write_zarr(adata, zarr_path, arr_path=["layers", "counts"], mode="w", client=client)

    return LazyAnnData(zarr_path, client=client)