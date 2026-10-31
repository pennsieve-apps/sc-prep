"""
Finding the samples under INPUT_DIR: Cell Ranger matrices (count, multi / Flex, ARC) by their
file and folder names, AnnData .h5ad and Seurat .rds files; the optional sample sheet; the SPARC
samples / subjects sheets. The files come from safe_inputs.InputDir.walk() (the platform's symlink
rule); every file read later goes through the same InputDir's checks or capped readers.
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from processor.safe_inputs import InputDir, UnsafeInput, assert_defusedxml

log = logging.getLogger("processor-sc-prep")

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
REQUIRED_COLUMNS = ("sample_id", "path", "source", "n_obs", "n_genes", "organism")
# Folders Cell Ranger names itself; a sample's default id is the first folder above them.
OUTS_NAMES = {"outs", "count"}
PORT = "counts"
SPARC_NAMES = {f"{stem}.{ext}" for stem in ("samples", "subjects") for ext in ("xlsx", "csv", "tsv")}
MTX_NAMES = ("matrix.mtx.gz", "matrix.mtx")
MAX_TABLE_BYTES = 50 * 1024 * 1024          # a sample sheet or SPARC table (CSV / TSV, or .xlsx on disk)
MAX_XLSX_UNCOMPRESSED = 200 * 1024 * 1024   # the declared size of an .xlsx's members
MAX_IGNORED = 500


def sanitise(stem: str) -> str:
    """A sample_id from a folder or file name: anything outside [A-Za-z0-9_-] becomes '_'."""
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_")
    return re.sub(r"_+", "_", s) or "sample"


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


@dataclass
class Found:
    candidates: list[Candidate] = field(default_factory=list)
    ignored: list[dict] = field(default_factory=list)
    n_ignored: int = 0
    sheet: Path | None = None
    sparc: dict[str, Path] = field(default_factory=dict)       # "samples" / "subjects" -> file
    n_files: int = 0


def default_id_for(folder: Path, root: Path) -> str:
    """The sample id without a sample sheet: the folder's name, or the first one above outs/ and count/
    (the input port's own folder, INPUT_DIR/counts/, does not name a sample)."""
    p = folder
    while p != root and (p.name in OUTS_NAMES or (p.name == PORT and p.parent == root)):
        p = p.parent
    return "sample" if p == root else sanitise(p.name)


def detect_folder(folder: Path, subdirs: list[str], names: list[str], is_mtx, indir: InputDir,
                  ignored: list[dict]) -> tuple[list[Candidate], bool]:
    """Candidates in one folder, and whether the folder is a Cell Ranger sample (its sub-folders are then not
    looked at). Raw matrices without a filtered one are added to `ignored`."""
    files = {n: folder / n for n in names}
    filtered = sorted(n for n in names if n.endswith("filtered_feature_bc_matrix.h5"))
    raw = sorted(n for n in names if n.endswith("raw_feature_bc_matrix.h5"))
    mtx = {d: folder / d for d in subdirs if d.endswith("feature_bc_matrix") and is_mtx(folder / d)}
    filtered += sorted(d for d in mtx if d.endswith("filtered_feature_bc_matrix"))
    raw += sorted(d for d in mtx if d.endswith("raw_feature_bc_matrix"))
    files.update(mtx)
    # an MTX folder given on its own (e.g. a dataset's aggregated filtered_feature_bc_matrix/ linked as the input)
    if not filtered and is_mtx(folder) and not folder.name.endswith("raw_feature_bc_matrix"):
        files[folder.name] = folder
        filtered = [folder.name]
    if filtered:
        c = Candidate("cellranger", folder, {"matrix": files[filtered[0]]}, default_id_for(folder, indir.root))
        if raw:
            c.files["raw"] = files[raw[0]]
        if "atac_fragments.tsv.gz" in names:
            c.files["atac_fragments"] = files["atac_fragments.tsv.gz"]
        return [c], True
    cands = []
    for n in raw:
        ignored.append({"path": indir.rel(files[n]), "reason": "raw matrix without a filtered matrix next to it"})
    for n in names:
        low = n.lower()
        if low.endswith(".h5ad"):
            cands.append(Candidate("h5ad", files[n], {"h5ad": files[n]}, sanitise(n[:-5])))
        elif low.endswith(".rds"):
            cands.append(Candidate("rds", files[n], {"rds": files[n]}, sanitise(n[:-4])))
    return cands, False


def walk(indir: InputDir, sheet_name: str) -> Found:
    """Every sample in the files InputDir.walk() accepts, the sample sheet and the SPARC sheets (shallowest
    wins; a file inside a Cell Ranger sample folder only when there is no other)."""
    found = Found()
    root = indir.root
    files = indir.walk()
    found.n_files = len(files)
    by_dir: dict[Path, list[str]] = {}
    for f in files:
        by_dir.setdefault(f.parent, []).append(f.name)
    children: dict[Path, set[str]] = {}
    for d in list(by_dir):
        while d != root:
            children.setdefault(d.parent, set()).add(d.name)
            d = d.parent
    is_mtx = lambda d: any(n in by_dir.get(d, ()) for n in MTX_NAMES)  # noqa: E731
    tables: dict[str, list[Path]] = {"sheet": [], "samples": [], "subjects": []}
    for f in files:
        if sheet_name and f.name == sheet_name:
            tables["sheet"].append(f)
        if f.name.lower() in SPARC_NAMES:
            tables[f.name.lower().split(".")[0]].append(f)
    claimed: list[Path] = []
    for folder in sorted(set(by_dir) | set(children), key=lambda p: (len(p.parts), str(p))):
        if any(folder.is_relative_to(c) for c in claimed):
            continue
        names = sorted(by_dir.get(folder, []))
        n_ignored = len(found.ignored)
        cands, whole = detect_folder(folder, sorted(children.get(folder, ())), names, is_mtx, indir, found.ignored)
        found.candidates.extend(cands)
        if whole:
            claimed.append(folder)
            continue
        used = {c.path for c in cands} | {root / i["path"] for i in found.ignored[n_ignored:]}
        for f in names:
            if folder / f not in used and f != sheet_name and f.lower() not in SPARC_NAMES:
                found.ignored.append({"path": indir.rel(folder / f), "reason": "not part of a recognised sample"})
    for key, paths in tables.items():
        paths = [p for p in paths if not any(p.is_relative_to(c) for c in claimed)] or paths
        if not paths:
            continue
        paths.sort(key=lambda p: (len(p.relative_to(root).parts), str(p)))
        if len(paths) > 1:
            log.warning("%d files for the %s table; using %s", len(paths), key, indir.rel(paths[0]))
        if key == "sheet":
            found.sheet = paths[0]
        else:
            found.sparc[key] = paths[0]
    found.candidates.sort(key=lambda c: indir.rel(c.path))
    found.ignored.sort(key=lambda i: i["path"])
    found.n_ignored = len(found.ignored)
    found.ignored = found.ignored[:MAX_IGNORED]
    return found


# --- tables -----------------------------------------------------------------------------------

def read_table(indir: InputDir, path: Path) -> pd.DataFrame:
    """A CSV, TSV or Excel sheet (first worksheet) as strings, empty cells as ''; header names stripped. Capped:
    UnsafeInput (a ValueError) when the file, or an .xlsx's uncompressed members, are too large."""
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        if os.path.getsize(indir.check(path)) > MAX_TABLE_BYTES:
            raise UnsafeInput(f"{indir.rel(path)} is larger than {MAX_TABLE_BYTES} bytes")
        if indir.xlsx_uncompressed_size(path) > MAX_XLSX_UNCOMPRESSED:
            raise UnsafeInput(f"{indir.rel(path)} is larger than {MAX_XLSX_UNCOMPRESSED} bytes uncompressed")
        assert_defusedxml()
        df = pd.read_excel(indir.check(path), sheet_name=0, dtype=str, engine="openpyxl")
    else:
        df = indir.read_csv(path, MAX_TABLE_BYTES, sep="\t" if suffix == ".tsv" else None, engine="python",
                            dtype=str, keep_default_na=False, encoding="utf-8-sig")
    df = df.fillna("")
    df.columns = [str(c).strip() for c in df.columns]
    df = df.loc[:, [not c.startswith("Unnamed:") for c in df.columns]]
    return df.apply(lambda s: s.str.strip())


@dataclass
class SampleSheet:
    path: Path | None = None
    columns: list[str] = field(default_factory=list)     # columns besides sample_id, path and the suffix column
    rows: list[dict[str, str]] = field(default_factory=list)


def read_sheet(indir: InputDir, path: Path | None, suffix_column: str) -> SampleSheet:
    """The sample sheet: needs sample_id and path columns. Raises ValueError when it cannot be used."""
    if path is None:
        return SampleSheet()
    df = read_table(indir, path)
    for col in ("sample_id", "path"):
        if col not in df.columns:
            raise ValueError(f"{path.name} has no '{col}' column (columns: {', '.join(df.columns)})")
    rows = [r for r in df.to_dict("records") if r["sample_id"] or r["path"]]
    extra = [c for c in df.columns if c not in REQUIRED_COLUMNS and c != suffix_column]
    return SampleSheet(path, extra, rows)


def _norm(p: str) -> str:
    return p.replace("\\", "/").strip().strip("/")


def path_matches(sheet_path: str, c: Candidate, indir: InputDir) -> bool:
    """A sheet `path` names a sample when it equals or ends the sample's path (with or without trailing outs/ or count/)."""
    want = _norm(sheet_path)
    if not want:
        return False
    forms = {_norm(indir.rel(c.path))}
    p = Path(indir.rel(c.path))
    while p.name in OUTS_NAMES and str(p.parent) != ".":
        p = p.parent
        forms.add(_norm(str(p)))
    return any(f == want or f.endswith("/" + want) for f in forms)


def assign_identities(cands: list[Candidate], sheet: SampleSheet, indir: InputDir, suffix_column: str) -> list[dict]:
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
        hits = [c for c in cands if path_matches(path, c, indir)]
        reason = None
        if not hits:
            reason = f"path {path!r} matches no sample under the input"
        elif len(hits) > 1:
            reason = f"path {path!r} matches {len(hits)} samples: " + ", ".join(indir.rel(h.path) for h in hits)
        elif suffix and hits[0].source != "cellranger":
            reason = f"{suffix_column} is set but {indir.rel(hits[0].path)} is not a Cell Ranger matrix"
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
