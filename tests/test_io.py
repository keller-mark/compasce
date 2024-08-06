import pytest

import zarr
import numpy as np
import pandas as pd
from os.path import join
from anndata import read_h5ad

import comparisce as csc

DATA_DIR = join("data")

def test_save():

    adata = read_h5ad(join(DATA_DIR, "lake_et_al.subset.h5ad"))
    zarr_path = join(DATA_DIR, "lake_et_al.subset.h5ad.zarr")

    client = csc.create_dask_client()

    ladata = csc.io.create_lazy_anndata(adata, zarr_path, client=client)
    ladata.save()

    # Check that things have been written properly.
    z = zarr.open(zarr_path, mode="r")

    assert str(z.tree()).strip() == """
/
 ├── X
 │   ├── data (3593680,) float32
 │   ├── indices (3593680,) int32
 │   └── indptr (20001,) int32
 ├── layers
 │   └── counts (20000, 3000) float32
 ├── obs
 │   ├── BMI
 │   │   ├── categories (4,) object
 │   │   └── codes (20000,) int8
 │   ├── ClusterClass
 │   │   ├── categories (4,) object
 │   │   └── codes (20000,) int8
 │   ├── ClusterNumber
 │   │   ├── categories (58,) int64
 │   │   └── codes (20000,) int8
 │   ├── SpecimenID
 │   │   ├── categories (47,) object
 │   │   └── codes (20000,) int8
 │   ├── assay
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── assay_ontology_term_id
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── cell_type
 │   │   ├── categories (26,) object
 │   │   └── codes (20000,) int8
 │   ├── cell_type_ontology_term_id
 │   │   ├── categories (26,) object
 │   │   └── codes (20000,) int8
 │   ├── celltype
 │   │   ├── categories (57,) object
 │   │   └── codes (20000,) int8
 │   ├── dataSource
 │   │   ├── categories (2,) object
 │   │   └── codes (20000,) int8
 │   ├── development_stage
 │   │   ├── categories (21,) object
 │   │   └── codes (20000,) int8
 │   ├── development_stage_ontology_term_id
 │   │   ├── categories (21,) object
 │   │   └── codes (20000,) int8
 │   ├── diabetes_history
 │   │   ├── categories (3,) object
 │   │   └── codes (20000,) int8
 │   ├── disease
 │   │   ├── categories (3,) object
 │   │   └── codes (20000,) int8
 │   ├── disease_ontology_term_id
 │   │   ├── categories (3,) object
 │   │   └── codes (20000,) int8
 │   ├── diseasetype
 │   │   ├── categories (3,) object
 │   │   └── codes (20000,) int8
 │   ├── donor_id
 │   │   ├── categories (45,) object
 │   │   └── codes (20000,) int8
 │   ├── eGFR
 │   │   ├── categories (5,) object
 │   │   └── codes (20000,) int8
 │   ├── hypertension
 │   │   ├── categories (2,) object
 │   │   └── codes (20000,) int8
 │   ├── index (20000,) object
 │   ├── is_primary_data (20000,) bool
 │   ├── nCount_RNA (20000,) float64
 │   ├── nFeature_RNA (20000,) int32
 │   ├── observation_joinid (20000,) object
 │   ├── organism
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── organism_ontology_term_id
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── percent.mt (20000,) float64
 │   ├── sampletype
 │   │   ├── categories (4,) object
 │   │   └── codes (20000,) int8
 │   ├── self_reported_ethnicity
 │   │   ├── categories (4,) object
 │   │   └── codes (20000,) int8
 │   ├── self_reported_ethnicity_ontology_term_id
 │   │   ├── categories (4,) object
 │   │   └── codes (20000,) int8
 │   ├── sex
 │   │   ├── categories (2,) object
 │   │   └── codes (20000,) int8
 │   ├── sex_ontology_term_id
 │   │   ├── categories (2,) object
 │   │   └── codes (20000,) int8
 │   ├── subclass.l1
 │   │   ├── categories (13,) object
 │   │   └── codes (20000,) int8
 │   ├── subclass.l2
 │   │   ├── categories (57,) object
 │   │   └── codes (20000,) int8
 │   ├── suspension_type
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── tissue
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── tissue_ontology_term_id
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   ├── tissue_type
 │   │   ├── categories (1,) object
 │   │   └── codes (20000,) int8
 │   └── tissuetype
 │       ├── categories (2,) object
 │       └── codes (20000,) int8
 ├── obsm
 │   ├── X_harmony (20000, 50) float64
 │   ├── X_pca (20000, 50) float64
 │   ├── X_umap1 (20000, 2) float64
 │   └── X_umap2 (20000, 2) float64
 ├── obsp
 │   └── distances
 │       ├── data (260968,) float64
 │       ├── indices (260968,) int32
 │       └── indptr (20001,) int32
 ├── uns
 │   ├── citation () <U302
 │   ├── default_embedding () <U7
 │   ├── schema_reference () <U87
 │   ├── schema_version () <U5
 │   ├── title () <U59
 │   └── write_metadata
 │       └── layers
 │           └── counts
 ├── var
 │   ├── _index (3000,) object
 │   ├── feature_biotype
 │   │   ├── categories (1,) object
 │   │   └── codes (3000,) int8
 │   ├── feature_is_filtered (3000,) bool
 │   ├── feature_length
 │   │   ├── categories (2543,) int64
 │   │   └── codes (3000,) int16
 │   ├── feature_name
 │   │   ├── categories (3000,) object
 │   │   └── codes (3000,) int16
 │   ├── feature_reference
 │   │   ├── categories (1,) object
 │   │   └── codes (3000,) int8
 │   ├── vst.mean (3000,) float64
 │   ├── vst.variable (3000,) int32
 │   ├── vst.variance (3000,) float64
 │   ├── vst.variance.expected (3000,) float64
 │   └── vst.variance.standardized (3000,) float64
 ├── varm
 │   ├── HARMONY (3000, 50) float64
 │   └── PCs (3000, 50) float64
 └── varp
 """.strip()
    