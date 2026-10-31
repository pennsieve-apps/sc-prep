"""
Readers. Each reads one input into raw AnnData pieces (counts as CSR cells x features, the
source's barcodes as obs_names, var with gene_ids / feature_types, any per-cell columns of the
source in obs) plus a `meta` dict for prep.json. processor/contract.py turns them into the
contract shape; nothing here filters or normalises. Also the droplets writer (raw matrices to
10x HDF5).
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from processor.inputs import Candidate, rel, sanitise

GEX = "Gene Expression"
SYMBOL_COLUMNS = ("feature_name", "gene_symbols", "gene_symbol", "gene_name", "gene_names", "symbol")
RDS_TIMEOUT_S = 1800
CHUNK_BARCODES = 50_000


@dataclass
class Read:
    adata: ad.AnnData
    meta: dict
    sub_id: str | None = None           # set when one input holds several samples
    warnings: list[str] = field(default_factory=list)


# --- 10x matrices -----------------------------------------------------------------------------

def read_10x_h5(path: Path) -> tuple[sp.csr_matrix, np.ndarray, pd.DataFrame, dict]:
    """A 10x HDF5 feature-barcode matrix (format 2+) as (CSR barcodes x features, barcodes, var, file attributes)."""
    with h5py.File(path, "r") as h:
        if "matrix" not in h:
            raise ValueError(f"{path.name}: not a 10x HDF5 matrix (no 'matrix' group)")
        m = h["matrix"]
        n_features, n_barcodes = (int(v) for v in m["shape"][:])
        # int32 indices read as such (not via an int64 copy): half the memory of the largest array
        idx_type = np.int32 if n_features < 2**31 else np.int64
        csc = sp.csc_matrix((m["data"][:], m["indices"].astype(idx_type)[:], m["indptr"][:]), shape=(n_features, n_barcodes))
        barcodes = m["barcodes"][:].astype(str)
        ft = m["features"]
        names = ft["name"][:].astype(str)
        var = pd.DataFrame({"gene_ids": ft["id"][:].astype(str) if "id" in ft else names,
                            "feature_types": ft["feature_type"][:].astype(str) if "feature_type" in ft else GEX},
                           index=pd.Index(names))
        if "genome" in ft:
            var["genome"] = ft["genome"][:].astype(str)
        attrs = {k: (v.decode() if isinstance(v, bytes) else str(v)) for k, v in h.attrs.items()
                 if k in ("chemistry_description", "software_version")}
    return csc.T.tocsr(copy=False), barcodes, var, attrs


def read_10x_mtx(folder: Path) -> tuple[sp.csr_matrix, np.ndarray, pd.DataFrame, dict]:
    """A 10x MTX folder (matrix.mtx[.gz], features.tsv[.gz] or genes.tsv[.gz], barcodes.tsv[.gz]), as read_10x_h5."""
    from scipy.io import mmread

    def pick(*names: str) -> Path:
        for n in names:
            if (folder / n).is_file():
                return folder / n
        raise ValueError(f"{folder.name}: no {names[0]}")

    mtx = pick("matrix.mtx.gz", "matrix.mtx")
    with (gzip.open if mtx.suffix == ".gz" else open)(mtx, "rb") as fh:
        X = sp.csr_matrix(mmread(fh, spmatrix=False).T)
    feats = pd.read_csv(pick("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"), sep="\t", header=None, dtype=str)
    barcodes = pd.read_csv(pick("barcodes.tsv.gz", "barcodes.tsv"), sep="\t", header=None, dtype=str)[0].to_numpy()
    names = feats[1].to_numpy() if feats.shape[1] > 1 else feats[0].to_numpy()
    var = pd.DataFrame({"gene_ids": feats[0].to_numpy(),
                        "feature_types": feats[2].to_numpy() if feats.shape[1] > 2 else GEX}, index=pd.Index(names))
    if X.shape != (len(barcodes), len(names)):
        raise ValueError(f"{folder.name}: matrix is {X.shape[::-1]} but {len(names)} features x {len(barcodes)} barcodes listed")
    return X, barcodes, var, {}


def read_10x(path: Path):
    return read_10x_mtx(path) if path.is_dir() else read_10x_h5(path)


def chemistry_of(description: str, feature_types: set[str]) -> str | None:
    d = description.lower()
    if "Peaks" in feature_types or "multiome" in d or "arc" in d.split():
        return "multiome"
    if "fixed rna" in d or "flex" in d:
        return "flex"
    if "5'" in d:
        return "5p"
    if "3'" in d:
        return "3p"
    return None


def read_cellranger(c: Candidate, root: Path, params: dict) -> list[Read]:
    """Cell Ranger count / multi (Flex) / ARC: the filtered matrix; an aggregated matrix is split by barcode suffix
    when the sample sheet maps suffixes to samples."""
    X, barcodes, var, attrs = read_10x(c.files["matrix"])
    types = set(var["feature_types"].astype(str))
    meta = {"source": "cellranger-arc" if "Peaks" in types else "cellranger",
            "files": {k: rel(v, root) for k, v in c.files.items()},
            "chemistry": chemistry_of(attrs.get("chemistry_description", ""), types),
            "software_version": attrs.get("software_version"),
            "genome": sorted(set(var["genome"])) if "genome" in var else None,
            "atac_present": "Peaks" in types, "atac_fragments": rel(c.files["atac_fragments"], root)
            if "atac_fragments" in c.files else None}
    var = var.drop(columns=["genome"], errors="ignore")
    if not c.suffix_rows:
        return [Read(ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(barcodes)), var=var), meta)]
    suffixes = np.array([b.rsplit("-", 1)[-1] if "-" in b else "" for b in barcodes])
    reads = []
    for suffix, row in c.suffix_rows.items():
        idx = np.flatnonzero(suffixes == suffix)
        reads.append(Read(ad.AnnData(X=X[idx], obs=pd.DataFrame(index=pd.Index(barcodes[idx])), var=var),
                          {**meta, "barcode_suffix": suffix}, sub_id=row["sample_id"] or sanitise(f"{c.default_id}_{suffix}")))
    unmapped = int(np.isin(suffixes, list(c.suffix_rows), invert=True).sum())
    if unmapped:
        reads[0].warnings.append(f"{unmapped} barcodes of {rel(c.files['matrix'], root)} have a suffix the sample "
                                 "sheet does not map; not used")
    return reads


# --- AnnData / Seurat -------------------------------------------------------------------------

def _integer_valued(M) -> bool:
    d = M.data if sp.issparse(M) else np.asarray(M)
    return bool(d.size == 0 or (np.all(np.isfinite(d)) and d.min() >= 0 and np.all(np.mod(d, 1) == 0)))


def _symbols(var: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """var with gene symbols as its index (from a symbol column when var_names are Ensembl ids), gene_ids, feature_types."""
    var = var.copy()
    names = var.index.astype(str)
    swapped = None
    if np.mean([n.startswith("ENS") for n in names[:1000]]) > 0.5:
        col = next((c for c in SYMBOL_COLUMNS if c in var.columns), None)
        if col is not None:
            var["gene_ids"] = names
            sym = var[col].astype(str).to_numpy()
            var.index = pd.Index(np.where((sym == "") | (sym == "nan"), names, sym))
            swapped = col
    if "gene_ids" not in var:
        var["gene_ids"] = np.where([n.startswith("ENS") for n in names], names, "")
    if "feature_types" not in var:
        var["feature_types"] = GEX
    return var[["gene_ids", "feature_types"]].astype(str), swapped


def read_h5ad_file(path: Path) -> tuple[ad.AnnData, str, list[str]]:
    """Counts from layers['counts'], else X, else raw.X (the first integer-valued one), with obs and var. Only these
    elements are read (not the whole file)."""
    warnings: list[str] = []
    with h5py.File(path, "r") as h:
        obs = ad.io.read_elem(h["obs"])
        var = ad.io.read_elem(h["var"])
        slots = [("layers/counts", "layers" in h and "counts" in h["layers"]), ("X", "X" in h),
                 ("raw/X", "raw" in h and "X" in h["raw"])]
        for slot, present in slots:
            if not present:
                continue
            M = ad.io.read_elem(h[slot])
            if _integer_valued(M):
                if slot == "raw/X":
                    var = ad.io.read_elem(h["raw/var"])
                break
            warnings.append(f"{slot} is not integer-valued (not raw counts)")
            del M
        else:
            raise ValueError("no raw counts: layers['counts'], X and raw.X are missing or not integer-valued")
    X = M.tocsr() if sp.issparse(M) else sp.csr_matrix(np.asarray(M))
    v, swapped = _symbols(var)
    if swapped:
        warnings.append(f"Ensembl var_names replaced by var['{swapped}'] symbols")
    obs.index = obs.index.astype(str)
    return ad.AnnData(X=X, obs=obs, var=v), slot, warnings


def _split(a: ad.AnnData, columns: tuple[str, ...]) -> tuple[str | None, list[tuple[str | None, ad.AnnData]]]:
    """One AnnData per value of the first of `columns` with more than one value (that column dropped), else as is."""
    for col in columns:
        if col in a.obs and a.obs[col].nunique(dropna=True) > 1:
            vals = a.obs[col].astype(str).to_numpy()
            parts = [(sanitise(v), a[vals == v].copy()) for v in pd.unique(vals)]
            for _, p in parts:
                p.obs = p.obs.drop(columns=[col])
            return col, parts
    return None, [(None, a)]


def _reads_from_anndata(a: ad.AnnData, meta: dict, warnings: list[str], split_columns: tuple[str, ...]) -> list[Read]:
    col, parts = _split(a, split_columns)
    if col:
        meta = {**meta, "split_by": col}
    return [Read(p, meta, sub_id=sid, warnings=list(warnings)) for sid, p in parts]


def read_h5ad(c: Candidate, root: Path, params: dict) -> list[Read]:
    a, slot, warnings = read_h5ad_file(c.files["h5ad"])
    meta = {"source": "h5ad", "files": {"h5ad": rel(c.files["h5ad"], root)}, "counts_from": slot}
    return _reads_from_anndata(a, meta, warnings, ("sample_id",))


def memory_limit_bytes() -> int:
    """80% of the container's memory (cgroup v2 / v1 limit, the Lambda function's size, else physical memory)."""
    for f in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            v = open(f).read().strip()
            if v.isdigit() and int(v) < 2**60:
                return int(int(v) * 0.8)
        except OSError:
            pass
    if os.environ.get("AWS_LAMBDA_FUNCTION_MEMORY_SIZE", "").isdigit():
        return int(int(os.environ["AWS_LAMBDA_FUNCTION_MEMORY_SIZE"]) * 2**20 * 0.8)
    return int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") * 0.8)


def read_rds(c: Candidate, root: Path, params: dict) -> list[Read]:
    """A Seurat object, converted by processor/rds.py in a child process (isolated, memory-capped, timed out)."""
    with tempfile.TemporaryDirectory(dir="/tmp" if os.path.isdir("/tmp") else None) as tmp:
        out = Path(tmp) / "counts.h5ad"
        env = {"PATH": "/usr/bin:/bin", "HOME": tmp, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
               "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1", "NUMBA_CACHE_DIR": tmp}
        try:
            proc = subprocess.run([sys.executable, "-I", "-W", "ignore", str(Path(__file__).with_name("rds.py")),
                                   str(c.files["rds"]), str(out), str(memory_limit_bytes())],
                                  env=env, cwd=tmp, capture_output=True, text=True, timeout=RDS_TIMEOUT_S)
        except subprocess.TimeoutExpired as exc:
            raise ValueError(f"reading the .rds took longer than {RDS_TIMEOUT_S} s") from exc
        if proc.returncode != 0:
            last = (proc.stderr.strip().splitlines() or [f"exit code {proc.returncode}"])[-1][:300]
            raise ValueError(f"cannot read the .rds: {last}")
        info = json.loads(proc.stdout.strip().splitlines()[-1])
        a, slot, warnings = read_h5ad_file(out)
    if info["assay"] != info["active_assay"]:
        warnings.append(f"read assay {info['assay']}; the object's active assay is {info['active_assay']}")
    meta = {"source": "rds", "files": {"rds": rel(c.files["rds"], root)}, "seurat": info}
    return _reads_from_anndata(a, meta, warnings, ("sample_id", "orig.ident"))


READERS = {"cellranger": read_cellranger, "h5ad": read_h5ad, "rds": read_rds}


# --- droplets ---------------------------------------------------------------------------------

def _write_10x_h5(dst: Path, barcodes: np.ndarray, var: pd.DataFrame, blocks, attrs: dict) -> int:
    """A 10x HDF5 matrix (features x barcodes, CSC) from blocks of (data, indices, per-barcode counts), in barcode
    order. Returns the number of non-zeros."""
    with h5py.File(dst, "w") as h:
        for k, v in {"filetype": "matrix", "version": 2, **attrs}.items():
            h.attrs[k] = v
        m = h.create_group("matrix")
        data = m.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(2**16,), compression="gzip")
        ind = m.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=(2**16,), compression="gzip")
        indptr = [np.zeros(1, dtype=np.int64)]
        n = 0
        for d, i, per in blocks:
            data.resize((n + len(d),))
            ind.resize((n + len(i),))
            data[n:] = d
            ind[n:] = i
            indptr.append(n + np.cumsum(per, dtype=np.int64))
            n += len(d)
        m.create_dataset("indptr", data=np.concatenate(indptr))
        m.create_dataset("shape", data=np.array([len(var), len(barcodes)], dtype=np.int32))
        m.create_dataset("barcodes", data=np.asarray(barcodes).astype("S"), compression="gzip")
        f = m.create_group("features")
        f.create_dataset("name", data=var.index.to_numpy().astype("S"))
        f.create_dataset("id", data=var["gene_ids"].to_numpy().astype("S"))
        f.create_dataset("feature_type", data=var["feature_types"].to_numpy().astype("S"))
        f.create_dataset("genome", data=(var["genome"] if "genome" in var else pd.Series([""] * len(var))).to_numpy().astype("S"))
        f.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
    return n


def write_droplets(raw: Path, dst: Path) -> dict:
    """The raw matrix as <sample_id>.raw.h5: copied when it is a 10x HDF5 of gene expression only; ARC's Peaks are
    left out (streamed by barcode blocks); an MTX folder is converted."""
    if raw.is_dir():
        X, barcodes, var, attrs = read_10x_mtx(raw)
        keep = np.flatnonzero(var["feature_types"].astype(str).to_numpy() == GEX)
        csc = X[:, keep].T.tocsc()
        n = _write_10x_h5(dst, barcodes, var.iloc[keep], [(csc.data.astype(np.int32), csc.indices, np.diff(csc.indptr))], attrs)
        return {"converted_from": "mtx", "n_barcodes": len(barcodes), "n_features": len(keep), "nnz": n,
                "features_dropped": int(len(var) - len(keep))}
    with h5py.File(raw, "r") as h:
        if "matrix" not in h:
            raise ValueError(f"{raw.name}: not a 10x HDF5 matrix")
        m = h["matrix"]
        types = m["features"]["feature_type"][:].astype(str) if "feature_type" in m["features"] else np.full(int(m["shape"][0]), GEX)
        n_barcodes = int(m["shape"][1])
        keep = types == GEX
        if keep.all():
            shutil.copyfile(raw, dst)
            return {"copied": True, "n_barcodes": n_barcodes, "n_features": int(keep.sum()), "features_dropped": 0}
        new_index = np.cumsum(keep) - 1
        ft = m["features"]
        var = pd.DataFrame({"gene_ids": ft["id"][:].astype(str)[keep], "feature_types": types[keep]},
                           index=pd.Index(ft["name"][:].astype(str)[keep]))
        if "genome" in ft:
            var["genome"] = ft["genome"][:].astype(str)[keep]
        indptr = m["indptr"][:]
        attrs = {k: v for k, v in h.attrs.items() if k in ("chemistry_description", "software_version")}

        def blocks():
            for s in range(0, n_barcodes, CHUNK_BARCODES):
                e = min(s + CHUNK_BARCODES, n_barcodes)
                lo, hi = int(indptr[s]), int(indptr[e])
                idx = m["indices"][lo:hi]
                d = m["data"][lo:hi]
                k = keep[idx]
                col = np.repeat(np.arange(e - s), np.diff(indptr[s:e + 1]))
                yield d[k].astype(np.int32), new_index[idx[k]], np.bincount(col[k], minlength=e - s)

        n = _write_10x_h5(dst, m["barcodes"][:], var, blocks(), attrs)
    return {"copied": False, "n_barcodes": n_barcodes, "n_features": int(keep.sum()), "nnz": n,
            "features_dropped": int((~keep).sum())}
