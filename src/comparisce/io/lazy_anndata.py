import zarr
import dask.array as da

from .zarr_io import dispatched_read_zarr, dispatched_write_zarr


class LazyAnnData:

    def __init__(self, zarr_path, var_chunk_size=5, client=None):
        self.zarr_path = zarr_path
        self.adata = dispatched_read_zarr(zarr_path)
        self.var_chunk_size = var_chunk_size
        self.client = client
    
    def __getitem__(self, key):
        print(f"Getting {key}")
        return self.adata.__dict__[key]
    
    def __setitem__(self, key, value):
        print(f"Setting {key}")
        self.adata.__dict__[key] = value
    
    def save(self, arr_path=None, mode="r+"):
        dispatched_write_zarr(self.adata, self.zarr_path, var_chunk_size=self.var_chunk_size, arr_path=arr_path, mode=mode)


    def get_da_from_zarr_layer(self, layer_key):
        X_path = f"/layers/{layer_key}"

        if self.client is None:
            raise ValueError("A Dask client must be provided to use get_da_from_zarr_layer")
            # TODO: instead of throwing, just run da.from_array() on the numpy array normally

        z = zarr.open(self.zarr_path, path=X_path, mode='r')
        X_shape = z.shape

        # Use Zarr so that X does not get loaded into memory
        # which would cause Dask to include it in the task graph,
        # which causes the task graph to be too large.
        X_dask = da.from_zarr(url=self.zarr_path, component=X_path, chunks=(X_shape[0], self.var_chunk_size), inline_array=False)
        return X_dask

    def put_da_to_zarr_layer(self, layer_key, X_dask):
        from dask.distributed import progress

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

    dispatched_write_zarr(adata, zarr_path, arr_path=["layers", "counts"], mode="w")

    return LazyAnnData(zarr_path, client=client)