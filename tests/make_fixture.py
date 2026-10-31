"""
A tiny synthetic input tree with one sample per supported layout, small enough for unit tests:

  python -m tests.make_fixture data/fixture [--layout linked|flat|both]

writes data/fixture/linked/counts/ (what a linked `counts` input looks like) and/or data/fixture/flat/
(the same files at the root of INPUT_DIR) with:

  samples.xlsx, subjects.xlsx   SPARC sheets: two samples of two subjects (a clashing column, two columns that
                                fall together in snake_case, an empty column, clinical columns)
  primary/sub-A/sam-A1/count/   Cell Ranger multi (Flex) per-sample outs: sample_filtered_feature_bc_matrix.h5 +
                                sample_raw_feature_bc_matrix.h5
  primary/sub-B/sam-B1/count/   the same, and a raw-only folder primary/sub-B/sam-B2/count/ (ignored)
  arc1/outs/                    Cell Ranger ARC: filtered + raw .h5 with Gene Expression and Peaks, atac_fragments
  run3/outs/                    Cell Ranger count (mouse) as MTX folders: filtered_feature_bc_matrix/ + raw_…/
  aggr/filtered_feature_bc_matrix/  an aggregated MTX with barcode suffixes -1 and -2 (one sample without a sheet)
  tab.h5ad                      counts in layers["counts"], log-normalised X, Ensembl var_names + feature_name
  multi.h5ad                    counts in raw.X only, obs["sample_id"] with two samples (split into m1, m2)
  seurat.rds                    a Seurat v5 object with split layers (tests/data, made by tests/data/make_rds.R)
  matrix_only.rds               a bare dgCMatrix (not a Seurat object: that sample fails, the run goes on)
  notes.txt                     a stray file (ignored)
"""
from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA = Path(__file__).resolve().parent / "data"
HUMAN = ["ACTB", "GAPDH", "SNAP25", "PRPH", "NEFH", "SCN10A", "RPL3", "RPS6", "MT-CO1", "MT-ND1", "TBCE", "TBCE"]
HUMAN_IDS = [f"ENSG{i:011d}" for i in range(1, len(HUMAN) + 1)]
MOUSE = ["Actb", "Gapdh", "Snap25", "Prph", "Nefh", "Rpl3", "mt-Co1", "Malat1"]
MOUSE_IDS = [f"ENSMUSG{i:011d}" for i in range(1, len(MOUSE) + 1)]
N_PEAKS = 5


def counts(n: int, g: int, seed: int, lam: float = 3.0) -> np.ndarray:
    return np.random.default_rng(seed).poisson(lam, size=(n, g)).astype(np.int32)


def barcodes(n: int, start: int = 0, suffix: str = "-1") -> list[str]:
    return [f"AAAC{i:08d}{suffix}" for i in range(start, start + n)]


def write_10x_h5(path: Path, X: np.ndarray, bcs: list[str], names: list[str], ids: list[str], types: list[str],
                 chemistry: str, genome: str = "GRCh38") -> None:
    """A 10x HDF5 matrix (format version 2): features x barcodes as CSC."""
    import h5py

    path.parent.mkdir(parents=True, exist_ok=True)
    M = sp.csc_matrix(X.T)
    with h5py.File(path, "w") as h:
        h.attrs.update({"filetype": "matrix", "version": 2, "software_version": "cellranger-8.0.1",
                        "chemistry_description": chemistry})
        m = h.create_group("matrix")
        m.create_dataset("data", data=M.data.astype(np.int32))
        m.create_dataset("indices", data=M.indices.astype(np.int64))
        m.create_dataset("indptr", data=M.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array(M.shape, dtype=np.int32))
        m.create_dataset("barcodes", data=np.array(bcs, dtype="S"))
        f = m.create_group("features")
        f.create_dataset("name", data=np.array(names, dtype="S"))
        f.create_dataset("id", data=np.array(ids, dtype="S"))
        f.create_dataset("feature_type", data=np.array(types, dtype="S"))
        f.create_dataset("genome", data=np.array([genome] * len(names), dtype="S"))
        f.create_dataset("_all_tag_keys", data=np.array([b"genome"]))


def write_mtx(folder: Path, X: np.ndarray, bcs: list[str], names: list[str], ids: list[str], types: list[str]) -> None:
    from scipy.io import mmwrite

    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / "matrix.mtx.gz", "wb") as fh:
        mmwrite(fh, sp.coo_matrix(X.T), field="integer")
    with gzip.open(folder / "features.tsv.gz", "wt") as fh:
        fh.writelines(f"{i}\t{n}\t{t}\n" for i, n, t in zip(ids, names, types))
    with gzip.open(folder / "barcodes.tsv.gz", "wt") as fh:
        fh.writelines(b + "\n" for b in bcs)


def flex_sample(folder: Path, n: int, seed: int) -> None:
    X = counts(n, len(HUMAN), seed)
    X[0] = 0                                       # one empty cell (MIN_COUNTS_PER_CELL drops it)
    gex = ["Gene Expression"] * len(HUMAN)
    write_10x_h5(folder / "sample_filtered_feature_bc_matrix.h5", X, barcodes(n), HUMAN, HUMAN_IDS, gex,
                 "Fixed RNA Profiling")
    raw = np.vstack([X, counts(3 * n, len(HUMAN), seed + 1, lam=0.2)])
    write_10x_h5(folder / "sample_raw_feature_bc_matrix.h5", raw, barcodes(4 * n), HUMAN, HUMAN_IDS, gex,
                 "Fixed RNA Profiling")


def sparc_sheets(root: Path) -> None:
    pd.DataFrame({
        "sample id": ["sam-A1", "sam-B1", "sam-B2"], "subject id": ["sub-A", "sub-B", "sub-B"],
        "sample experimental group": ["", "", ""], "sample anatomical location": ["dorsal root ganglion"] * 3,
        "laterality": ["Left", "Right", "Right"], "Freezing method": ["dry ice", "", ""],
        "freezing method": ["", "liquid nitrogen", ""], "laboratory internal id": ["L1", "L2", "L3"],
    }).to_excel(root / "samples.xlsx", index=False)
    pd.DataFrame({
        "subject id": ["sub-A", "sub-B"], "subject experimental group": ["acute", "chronic"],
        "age": ["61 years", "70 years"], "sex": ["male", "female"], "species": ["Homo sapiens"] * 2,
        "laboratory internal id": ["P1", "P2"], "left_average_pain": ["7", "2"], "metadata only": ["", ""],
    }).to_excel(root / "subjects.xlsx", index=False)


def make_tree(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    sparc_sheets(root)
    flex_sample(root / "primary/sub-A/sam-A1/count", 30, 1)
    flex_sample(root / "primary/sub-B/sam-B1/count", 24, 2)
    b2 = root / "primary/sub-B/sam-B2/count"
    b2.mkdir(parents=True)
    shutil.copy(root / "primary/sub-B/sam-B1/count/sample_raw_feature_bc_matrix.h5", b2)

    # ARC: Gene Expression + Peaks in both matrices
    names = HUMAN + [f"chr1:{i * 1000}-{i * 1000 + 500}" for i in range(N_PEAKS)]
    ids = HUMAN_IDS + names[len(HUMAN):]
    types = ["Gene Expression"] * len(HUMAN) + ["Peaks"] * N_PEAKS
    X = counts(20, len(names), 3)
    outs = root / "arc1/outs"
    write_10x_h5(outs / "filtered_feature_bc_matrix.h5", X, barcodes(20), names, ids, types,
                 "Single Cell Multiome ATAC + Gene Expression v1")
    write_10x_h5(outs / "raw_feature_bc_matrix.h5", np.vstack([X, counts(40, len(names), 4, lam=0.2)]), barcodes(60),
                 names, ids, types, "Single Cell Multiome ATAC + Gene Expression v1")
    with gzip.open(outs / "atac_fragments.tsv.gz", "wt") as fh:
        fh.write("chr1\t100\t200\tAAAC00000000-1\t1\n")

    # Cell Ranger count, mouse, MTX folders
    gex = ["Gene Expression"] * len(MOUSE)
    X = counts(18, len(MOUSE), 5)
    write_mtx(root / "run3/outs/filtered_feature_bc_matrix", X, barcodes(18), MOUSE, MOUSE_IDS, gex)
    write_mtx(root / "run3/outs/raw_feature_bc_matrix", np.vstack([X, counts(30, len(MOUSE), 6, lam=0.2)]),
              barcodes(48), MOUSE, MOUSE_IDS, gex)

    # aggregated MTX, two libraries by barcode suffix
    gex = ["Gene Expression"] * len(HUMAN)
    write_mtx(root / "aggr/filtered_feature_bc_matrix", counts(22, len(HUMAN), 7),
              barcodes(12, suffix="-1") + barcodes(10, suffix="-2"), HUMAN, HUMAN_IDS, gex)

    import anndata as ad

    C = counts(16, len(HUMAN), 8)
    a = ad.AnnData(X=sp.csr_matrix(np.log1p(C).astype(np.float32)),
                   obs=pd.DataFrame({"cell_type": ["neuron", "glia"] * 8}, index=barcodes(16)),
                   var=pd.DataFrame({"feature_name": HUMAN}, index=HUMAN_IDS))
    a.layers["counts"] = sp.csr_matrix(C)
    a.write_h5ad(root / "tab.h5ad")
    C = counts(14, len(HUMAN[:10]), 9)
    full = ad.AnnData(X=sp.csr_matrix(C), obs=pd.DataFrame({"sample_id": ["m1"] * 7 + ["m2"] * 7}, index=barcodes(14)),
                      var=pd.DataFrame(index=HUMAN[:10]))
    m = full.copy()
    m.raw = full
    m.X = sp.csr_matrix(np.log1p(C).astype(np.float32))
    m.write_h5ad(root / "multi.h5ad")

    shutil.copy(DATA / "seurat_v5_split.rds", root / "seurat.rds")
    shutil.copy(DATA / "dgc.rds", root / "matrix_only.rds")
    (root / "notes.txt").write_text("not data\n")


def make_fixture(out: Path, layout: str = "both") -> dict[str, Path]:
    out = Path(out)
    roots = {}
    if layout in ("linked", "both"):
        make_tree(out / "linked" / "counts")
        roots["linked"] = out / "linked"
    if layout in ("flat", "both"):
        make_tree(out / "flat")
        roots["flat"] = out / "flat"
    return roots


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    args = ap.parse_args()
    for layout, path in make_fixture(args.out, args.layout).items():
        print(f"{layout}: {path}")
