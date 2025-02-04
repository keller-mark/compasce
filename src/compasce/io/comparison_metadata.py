import json
import numpy as np
import zarr
from os.path import join
from .cdata import dir_name_to_str


# Reference: https://stackoverflow.com/a/57915246
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


# See https://observablehq.com/d/e7c03bf319f20f86
class ComparisonMetadata:
    def __init__(self, comparison_key):
        self.comparison_key = comparison_key
        self.comparison_key_str = dir_name_to_str(comparison_key)
        self.items = []
    
    def get_df_key(self, df_type):
        return f"{self.comparison_key_str}.{df_type}"
    
    def append_df(self, adata_key, df_type, df_params, df_c_vals):
        self.items.append({
            "path": join(adata_key, self.get_df_key(df_type)),
            "coordination_values": df_c_vals,
            "analysis_type": df_type,
            "analysis_params": df_params,
        })
        return self.get_df_key(df_type)
    
    def get_dict(self):
        return {
            self.comparison_key_str: {
                "comparison": self.comparison_key,
                "results": self.items,
            }
        }
    
class MultiComparisonMetadata:
    def __init__(self, sample_group_pairs=None, sample_id_col=None, cell_type_col=None):
        self.schema_version = "0.0.1"
        self._comparisons = []
        self.sample_group_pairs = sample_group_pairs
        self.sample_id_col = sample_id_col
        self.cell_type_col = cell_type_col
        self._prev_comparisons_dict = dict()
    
    def load_state(self, zarr_path):
        z = zarr.open(zarr_path, mode="r+")
        if "/uns/comparison_metadata" in z:
            prev = json.loads(str(z["/uns/comparison_metadata"][()]))
            self._prev_comparisons_dict = prev["comparisons"]

    def add_comparison(self, comparison_key):
        c = ComparisonMetadata(comparison_key)
        self._comparisons.append(c)
        return c
    
    def serialize(self):
        comparisons_dict = self._prev_comparisons_dict
        for c in self._comparisons:
            comparisons_dict.update(c.get_dict())
        return json.dumps({
            "schema_version": self.schema_version,
            "comparisons": comparisons_dict,

            "sample_id_col": self.sample_id_col,
            "sample_group_pairs": self.sample_group_pairs,
            "cell_type_col": self.cell_type_col,
        })
    