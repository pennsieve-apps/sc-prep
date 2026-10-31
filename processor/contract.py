"""
From a reader's raw AnnData to the contract (docs/conventions.md in workflow-single-cell, "The
AnnData contract"): Gene Expression features only, gene symbols as unique var_names with mt /
ribo flags, integer counts as CSR, the optional count floor and seeded subsample, QC metrics,
obs_names prefixed with the sample id, sample metadata as categorical obs columns, and
uns["pennsieve"].
"""
from __future__ import annotations

import re
import zlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from processor.readers import GEX, Read

QC_OBS = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
# Per-cell columns of a source .h5ad / Seurat object that are labels (author cell types, clusters, the sample
# column): copied onto cells. Any other source column (author objects often carry clinical fields per cell) is
# copied only when named in SOURCE_COLUMNS; the names (never the values) are listed in source_columns.csv.
LABEL_RE = re.compile(r"^(orig_ident|class|.*_snn_res_.*|(.*_)?(leiden|louvain|clusters?|cell_?types?|annotations?|"
                      r"sub_?class(es)?|sub_?types?|labels?)(_.*)?)$", re.IGNORECASE)
# obs columns this and later steps write; a source column of the same name is kept as source_<name>
RESERVED = {"sample_id", "condition", *QC_OBS, "doublet_score", "predicted_doublet", "ambient_fraction", "cluster",
            "cell_type", "cell_type_score", "cell_type_source"}
MT_RE = re.compile(r"^mt-", re.IGNORECASE)
RIBO_RE = re.compile(r"^rp[sl]", re.IGNORECASE)
SPECIES = {"homo sapiens": "human", "human": "human", "mus musculus": "mouse", "mouse": "mouse",
           "rattus norvegicus": "rat", "rat": "rat"}
GENOMES = (("grch", "human"), ("hg", "human"), ("grcm", "mouse"), ("mm", "mouse"), ("rnor", "rat"),
           ("mratbn", "rat"), ("grcr", "rat"))
ENSEMBL = (("ENSMUSG", "mouse"), ("ENSRNOG", "rat"), ("ENSG", "human"))


def as_counts(X) -> sp.csr_matrix:
    """CSR int32 (int32 indices when they fit). Raises ValueError when the values are not non-negative integers."""
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    X.sum_duplicates()
    d = X.data
    if d.size and (not np.all(np.isfinite(d)) or d.min() < 0 or np.any(np.mod(d, 1) != 0)):
        raise ValueError("the matrix is not raw counts (negative or non-integer values)")
    if d.size and d.max() >= 2**31:
        raise ValueError("counts above 2^31")
    X = X.astype(np.int32, copy=False)
    if X.nnz < 2**31 - 1:
        X.indices = X.indices.astype(np.int32, copy=False)
        X.indptr = X.indptr.astype(np.int32, copy=False)
    return X


def subsample_index(n: int, max_obs: int, seed: int, sample_id: str) -> np.ndarray | None:
    """A seeded subset of row indices in source order (each sample its own subset); None when nothing is dropped."""
    if max_obs <= 0 or n <= max_obs:
        return None
    rng = np.random.default_rng([seed, zlib.crc32(sample_id.encode())])
    return np.sort(rng.choice(n, size=max_obs, replace=False))


def organism_of(param: str, species: str, genomes: list[str] | None, var: pd.DataFrame) -> tuple[str, str]:
    """(organism, how it was decided): the ORGANISM parameter, the metadata's species, the 10x genome, the
    Ensembl id prefix, else gene symbol case (all upper case: human)."""
    if param != "auto":
        return param, "parameter"
    if species.strip().lower() in SPECIES:
        return SPECIES[species.strip().lower()], "metadata species"
    for g in genomes or []:
        for prefix, org in GENOMES:
            if g.lower().startswith(prefix):
                return org, "genome"
    ids = var["gene_ids"].astype(str)
    for prefix, org in ENSEMBL:
        if ids.str.startswith(prefix).mean() > 0.5:
            return org, "Ensembl ids"
    names = [n for n in var.index.astype(str)[:2000] if n[:1].isalpha()]
    if names and np.mean([n.isupper() for n in names]) > 0.8:
        return "human", "gene symbol case"
    return "unknown", "undetermined (set ORGANISM)"


def finalise(read: Read, sample_id: str, source_rel: str, params: dict, metadata: dict[str, str],
             obs_columns: list[str], species: str, prep_version: str) -> tuple[ad.AnnData, dict]:
    """The contract AnnData and the prep.json record for one sample. `obs_columns` are the metadata columns put on
    every cell (values from `metadata`, '' as missing)."""
    a, meta = read.adata, read.meta
    warnings = list(read.warnings)
    var = a.var.copy()
    var.index = var.index.astype(str)
    n_obs_source, n_feat_source = a.shape

    # 1. features: Gene Expression only
    types = var["feature_types"].astype(str).to_numpy()
    gene = types == GEX
    if not gene.any():
        raise ValueError("no Gene Expression features")
    dropped = {t: int((types == t).sum()) for t in sorted(set(types[~gene]))}
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(np.asarray(a.X))
    if not gene.all():
        X = X[:, np.flatnonzero(gene)]
        var = var.iloc[np.flatnonzero(gene)]

    # 2. cells: count floor, then the subsample (before any copy of the whole matrix)
    obs = a.obs.copy()
    totals = np.asarray(X.sum(axis=1)).ravel()
    idx = np.flatnonzero(totals >= params["MIN_COUNTS_PER_CELL"])
    n_after_floor = len(idx)
    sub = subsample_index(n_after_floor, params["MAX_CELLS_PER_SAMPLE"], params["SEED"], sample_id)
    if sub is not None:
        idx = idx[sub]
    if len(idx) == 0:
        raise ValueError(f"no cells left (MIN_COUNTS_PER_CELL={params['MIN_COUNTS_PER_CELL']})")
    if len(idx) < n_obs_source:
        X, obs = X[idx], obs.iloc[idx]
    del a
    read.adata = None
    X = as_counts(X)

    # 3. var: unique symbols, mt / ribo
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=obs.index.astype(str)), var=var[["gene_ids", "feature_types"]].copy())
    n_dup = int(adata.var_names.duplicated().sum())
    if n_dup:
        adata.var_names_make_unique()
        warnings.append(f"{n_dup} duplicated gene symbols made unique")
    names = adata.var_names.to_series()
    adata.var["mt"] = names.str.match(MT_RE).to_numpy()
    adata.var["ribo"] = names.str.match(RIBO_RE).to_numpy()

    # 4. QC metrics
    import scanpy as sc

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
    adata.var = adata.var[["gene_ids", "feature_types", "mt", "ribo"]]
    qc = adata.obs[QC_OBS].copy()
    qc["total_counts"] = qc["total_counts"].astype(np.float64)
    qc["n_genes_by_counts"] = qc["n_genes_by_counts"].astype(np.int64)

    # 5. obs: ids, QC, metadata, then the source's own columns
    pre = f"{sample_id}:"
    barcodes = [o[len(pre):] if o.startswith(pre) else o for o in adata.obs_names]
    new = pd.DataFrame(index=pd.Index([pre + b for b in barcodes]))
    new["sample_id"] = pd.Categorical([sample_id] * adata.n_obs)
    for col in QC_OBS:
        new[col] = qc[col].to_numpy()
    for col in obs_columns:
        val = metadata.get(col, "")
        new[col] = pd.Categorical([val] * adata.n_obs if val else [np.nan] * adata.n_obs)
    renamed, source_columns = [], {}
    for col in obs.columns:
        if not (LABEL_RE.match(re.sub(r"[.\-\s]+", "_", str(col))) or str(col) in params["SOURCE_COLUMNS"]):
            source_columns[str(col)] = None
            continue
        name = f"source_{col}" if col in RESERVED or col in new.columns else str(col)
        source_columns[str(col)] = name
        if name != col:
            renamed.append(str(col))
        v = obs[col].to_numpy()
        new[name] = pd.Categorical(pd.Series(v).astype(str).to_numpy()) if obs[col].dtype == object or \
            isinstance(obs[col].dtype, (pd.CategoricalDtype, pd.StringDtype)) or str(obs[col].dtype) == "str" else v
    if renamed:
        warnings.append(f"source obs columns renamed to source_<name>: {', '.join(renamed)}")
    adata.obs = new

    organism, how = organism_of(params["ORGANISM"], species, meta.get("genome"), adata.var)
    if organism == "unknown":
        warnings.append("organism could not be determined; set ORGANISM")
    pen = {"source": source_rel, "source_type": meta["source"], "organism": organism, "prep_version": prep_version,
           "sample_id": sample_id, "n_obs_source": int(n_obs_source)}
    if meta.get("chemistry"):
        pen["chemistry"] = meta["chemistry"]
    if sub is not None:
        pen["subsampled_from"] = int(n_after_floor)
    adata.uns["pennsieve"] = pen

    info = {
        "source": meta["source"], "files": meta.get("files", {}), "chemistry": meta.get("chemistry"),
        "organism": organism, "organism_from": how,
        "features": {"n_source": int(n_feat_source), "n_genes": int(adata.n_vars), "dropped": dropped},
        "cells": {"n_source": int(n_obs_source), "n_below_min_counts": int(n_obs_source - n_after_floor),
                  "n_subsampled_to": int(len(idx)) if sub is not None else None, "n_written": int(adata.n_obs)},
        "atac_present": bool(meta.get("atac_present")), "n_obs": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "source_columns": source_columns,
        "warnings": warnings,
    }
    for k in ("software_version", "atac_fragments", "counts_from", "split_by", "barcode_suffix", "seurat"):
        if meta.get(k) is not None:
            info[k] = meta[k]
    return adata, info
