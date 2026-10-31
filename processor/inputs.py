"""
Finding the samples under INPUT_DIR: Cell Ranger matrices (count, multi / Flex, ARC) by their
file and folder names, AnnData .h5ad and Seurat .rds files; the optional sample sheet; the SPARC
samples / subjects sheets. Paths that resolve outside the input are refused.
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

log = logging.getLogger("processor-sc-prep")

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
REQUIRED_COLUMNS = ("sample_id", "path", "source", "n_obs", "n_genes", "organism")
# Folders Cell Ranger names itself; a sample's default id is the first folder above them.
OUTS_NAMES = {"outs", "count"}
PORT = "counts"
SPARC_NAMES = {f"{stem}.{ext}" for stem in ("samples", "subjects") for ext in ("xlsx", "csv", "tsv")}
MAX_TABLE_BYTES = 50 * 1024 * 1024
MAX_IGNORED = 500


def sanitise(stem: str) -> str:
    """A sample_id from a folder or file name: anything outside [A-Za-z0-9_-] becomes '_'."""
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_")
    return re.sub(r"_+", "_", s) or "sample"


def rel(path: Path, root: Path) -> str:
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


@dataclass
class Candidate:
    """One detected input: a Cell Ranger sample folder, or an .h5ad / .rds file."""
    source: str                        # cellranger | h5ad | rds (cellranger-arc is told apart when read)
    path: Path                         # the sample's folder, or its file
    files: dict[str, Path] = field(default_factory=dict)
    default_id: str = ""
    notes: list[str] = field(default_factory=list)
    sample_id: str = ""
    sheet_row: dict[str, str] | None = None
    suffix_rows: dict[str, dict[str, str]] = field(default_factory=dict)   # barcode suffix -> sheet row

    def rel(self, root: Path) -> str:
        return rel(self.path, root)


@dataclass
class Found:
    candidates: list[Candidate] = field(default_factory=list)
    ignored: list[dict] = field(default_factory=list)
    n_ignored: int = 0
    sheet: Path | None = None
    sparc: dict[str, Path] = field(default_factory=dict)       # "samples" / "subjects" -> file
    n_files: int = 0


def allowed_roots(root: Path) -> list[str]:
    """
    The real paths files may resolve to: INPUT_DIR itself, plus the targets of the symlinks the
    platform places in it (at its root, or in a port's folder), which point at earlier steps'
    outputs and the run's staged data. A symlink deeper in the data must stay inside these.
    """
    roots = [os.path.realpath(root)]
    for entry in root.iterdir():
        if entry.is_symlink():
            roots.append(os.path.realpath(entry))
        elif entry.is_dir():
            roots.extend(os.path.realpath(e) for e in entry.iterdir() if e.is_symlink())
    return roots


def inside(path: Path, roots: list[str]) -> bool:
    real = os.path.realpath(path)
    return any(real == r or real.startswith(r.rstrip(os.sep) + os.sep) for r in roots)


def default_id_for(folder: Path, root: Path) -> str:
    """The sample id without a sample sheet: the folder's name, or the first one above outs/ and count/
    (the input port's own folder, INPUT_DIR/counts/, does not name a sample)."""
    p = folder
    while p != root and (p.name in OUTS_NAMES or (p.name == PORT and p.parent == root)):
        p = p.parent
    return "sample" if p == root else sanitise(p.name)


def _is_mtx_dir(d: Path) -> bool:
    return d.is_dir() and any((d / n).is_file() for n in ("matrix.mtx.gz", "matrix.mtx"))


def detect_folder(folder: Path, dirnames: list[str], filenames: list[str], root: Path,
                  ignored: list[dict]) -> tuple[list[Candidate], bool]:
    """Candidates in one folder, and whether the folder is a Cell Ranger sample (not descended into). Raw matrices
    without a filtered one are added to `ignored`."""
    files = {n: folder / n for n in filenames}
    filtered = sorted(n for n in filenames if n.endswith("filtered_feature_bc_matrix.h5"))
    raw = sorted(n for n in filenames if n.endswith("raw_feature_bc_matrix.h5"))
    mtx = {d: folder / d for d in dirnames if d.endswith("feature_bc_matrix") and _is_mtx_dir(folder / d)}
    filtered += sorted(d for d in mtx if d.endswith("filtered_feature_bc_matrix"))
    raw += sorted(d for d in mtx if d.endswith("raw_feature_bc_matrix"))
    files.update(mtx)
    # an MTX folder given on its own (e.g. a dataset's aggregated filtered_feature_bc_matrix/ linked as the input)
    if not filtered and _is_mtx_dir(folder) and not folder.name.endswith("raw_feature_bc_matrix"):
        files[folder.name] = folder
        filtered = [folder.name]
    if filtered:
        c = Candidate("cellranger", folder, {"matrix": files[filtered[0]]}, default_id_for(folder, root))
        if raw:
            c.files["raw"] = files[raw[0]]
        if "atac_fragments.tsv.gz" in filenames:
            c.files["atac_fragments"] = files["atac_fragments.tsv.gz"]
        return [c], True
    cands = []
    for n in raw:
        ignored.append({"path": rel(files[n], root), "reason": "raw matrix without a filtered matrix next to it"})
    for n in filenames:
        low = n.lower()
        if low.endswith(".h5ad"):
            cands.append(Candidate("h5ad", files[n], {"h5ad": files[n]}, sanitise(n[:-5])))
        elif low.endswith(".rds"):
            cands.append(Candidate("rds", files[n], {"rds": files[n]}, sanitise(n[:-4])))
    return cands, False


def walk(root: Path, sheet_name: str) -> Found:
    """Every sample under root (symlinks followed inside the input), the sample sheet and SPARC sheets (shallowest wins)."""
    found = Found()
    roots = allowed_roots(root)
    seen: set[str] = set()
    tables: dict[str, list[Path]] = {"sheet": [], "samples": [], "subjects": []}
    claimed: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        folder = Path(dirpath)
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        if not inside(folder, roots):
            found.ignored.append({"path": rel(folder, root) + "/", "reason": "resolves outside the input; not read"})
            dirnames[:] = []
            continue
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        names = []
        for f in sorted(filenames):
            p = folder / f
            if f.startswith(".") or not p.is_file():
                continue
            if not inside(p, roots):
                found.ignored.append({"path": rel(p, root), "reason": "resolves outside the input; not read"})
                continue
            names.append(f)
        found.n_files += len(names)
        for f in names:
            if sheet_name and f == sheet_name:
                tables["sheet"].append(folder / f)
            if f.lower() in SPARC_NAMES:
                tables[f.lower().split(".")[0]].append(folder / f)
        n_ignored = len(found.ignored)
        cands, whole = detect_folder(folder, dirnames, names, root, found.ignored)
        found.candidates.extend(cands)
        if whole:
            claimed.append(folder)
            dirnames[:] = []
            continue
        used = {c.path for c in cands} | {root / i["path"] for i in found.ignored[n_ignored:]}
        for f in names:
            if folder / f not in used and f != sheet_name and f.lower() not in SPARC_NAMES:
                found.ignored.append({"path": rel(folder / f, root), "reason": "not part of a recognised sample"})
    for key, paths in tables.items():
        paths = [p for p in paths if not any(p.is_relative_to(c) for c in claimed)] or paths
        if not paths:
            continue
        paths.sort(key=lambda p: (len(p.relative_to(root).parts), str(p)))
        if len(paths) > 1:
            log.warning("%d files for the %s table; using %s", len(paths), key, rel(paths[0], root))
        if key == "sheet":
            found.sheet = paths[0]
        else:
            found.sparc[key] = paths[0]
    # a SPARC samples.csv with no `path` column is SPARC metadata, not the sample sheet
    if found.sheet is not None and found.sheet == found.sparc.get("samples"):
        if "path" in read_table(found.sheet).columns:
            del found.sparc["samples"]
        else:
            found.sheet = None
    found.candidates.sort(key=lambda c: c.rel(root))
    found.ignored.sort(key=lambda i: i["path"])
    found.n_ignored = len(found.ignored)
    found.ignored = found.ignored[:MAX_IGNORED]
    return found


# --- tables -----------------------------------------------------------------------------------

def read_table(path: Path) -> pd.DataFrame:
    """A CSV, TSV or Excel sheet (first worksheet) as strings, empty cells as ''. Header names stripped."""
    if path.stat().st_size > MAX_TABLE_BYTES:
        raise ValueError(f"{path.name} is larger than {MAX_TABLE_BYTES // 2**20} MB")
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        df = pd.read_excel(path, sheet_name=0, dtype=str, engine="openpyxl")
    else:
        df = pd.read_csv(path, sep="\t" if suffix == ".tsv" else None, engine="python", dtype=str,
                         keep_default_na=False, encoding="utf-8-sig")
    df = df.fillna("")
    df.columns = [str(c).strip() for c in df.columns]
    df = df.loc[:, [not c.startswith("Unnamed:") for c in df.columns]]
    return df.apply(lambda s: s.str.strip())


@dataclass
class SampleSheet:
    path: Path | None = None
    columns: list[str] = field(default_factory=list)     # columns besides sample_id, path and the suffix column
    rows: list[dict[str, str]] = field(default_factory=list)


def read_sheet(path: Path | None, suffix_column: str) -> SampleSheet:
    """The sample sheet: needs sample_id and path columns. Raises ValueError when it cannot be used."""
    if path is None:
        return SampleSheet()
    df = read_table(path)
    for col in ("sample_id", "path"):
        if col not in df.columns:
            raise ValueError(f"{path.name} has no '{col}' column (columns: {', '.join(df.columns)})")
    rows = [r for r in df.to_dict("records") if r["sample_id"] or r["path"]]
    extra = [c for c in df.columns if c not in REQUIRED_COLUMNS and c != suffix_column]
    return SampleSheet(path, extra, rows)


def _norm(p: str) -> str:
    return p.replace("\\", "/").strip().strip("/")


def path_matches(sheet_path: str, c: Candidate, root: Path) -> bool:
    """A sheet `path` names a sample when it equals or ends the sample's path (with or without trailing outs/ or count/)."""
    want = _norm(sheet_path)
    if not want:
        return False
    forms = {_norm(c.rel(root))}
    p = Path(c.rel(root))
    while p.name in OUTS_NAMES and str(p.parent) != ".":
        p = p.parent
        forms.add(_norm(str(p)))
    return any(f == want or f.endswith("/" + want) for f in forms)


def assign_identities(cands: list[Candidate], sheet: SampleSheet, root: Path, suffix_column: str) -> list[dict]:
    """
    Give each candidate its sample_id and sheet row. A row with a barcode suffix (suffix_column) names one
    sample inside an aggregated Cell Ranger matrix; several such rows may name the same matrix.
    Returns the rows that matched nothing, several samples, or a sample already taken, as error records.
    """
    errors: list[dict] = []
    for c in cands:
        c.sample_id = c.default_id
    for row in sheet.rows:
        sid, path, suffix = row.get("sample_id", ""), row.get("path", ""), row.get(suffix_column, "")
        hits = [c for c in cands if path_matches(path, c, root)]
        reason = None
        if not hits:
            reason = f"path {path!r} matches no sample under the input"
        elif len(hits) > 1:
            reason = f"path {path!r} matches {len(hits)} samples: " + ", ".join(h.rel(root) for h in hits)
        elif suffix and hits[0].source != "cellranger":
            reason = f"{suffix_column} is set but {hits[0].rel(root)} is not a Cell Ranger matrix"
        elif suffix and (hits[0].sheet_row is not None or suffix in hits[0].suffix_rows):
            reason = f"path {path!r} with {suffix_column} {suffix!r} is listed twice"
        elif not suffix and (hits[0].sheet_row is not None or hits[0].suffix_rows):
            reason = f"path {path!r} matches the same sample as an earlier row"
        if reason:
            errors.append({"sample_id": sid, "path": path, "reason": reason})
            continue
        c = hits[0]
        if sid and not SAMPLE_ID_RE.match(sid):
            c.notes.append(f"sample_id {sid!r} sanitised to {sanitise(sid)!r}")
            row = {**row, "sample_id": sanitise(sid)}
        if suffix:
            c.suffix_rows[suffix] = row
        else:
            c.sheet_row = row
            c.sample_id = row["sample_id"] or c.default_id
    return errors
