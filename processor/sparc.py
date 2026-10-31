"""
SPARC Dataset Structure metadata: samples.xlsx / subjects.xlsx (or .csv / .tsv) joined onto the
samples found. A sample folder named `sam-…` (or below one) matches the `sample id` row, else the
sample_id itself does; its `subject id` (or a `sub-…` folder above it) joins the subject row.
Column names become snake_case; subject columns get a `subject_` prefix only when they clash with
a sample column; columns empty in every row are dropped. Values are never logged.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from processor.inputs import read_table
from processor.safe_inputs import InputDir

KEYS = ("sample_id", "subject_id")


def snake(name: str) -> str:
    """`subject experimental group` -> `subject_experimental_group`, `wasDerivedFrom` -> `was_derived_from`."""
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", str(name).strip())
    s = re.sub(r"[^0-9a-zA-Z]+", "_", s).strip("_").lower()
    return s or "column"


def _snake_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Columns renamed to snake_case; columns that fall together (`Freezing method`, `freezing method`) are
    coalesced, first non-empty value wins; columns empty in every row dropped."""
    out: dict[str, pd.Series] = {}
    for col in df.columns:
        key = snake(col)
        vals = df[col].astype(str)
        out[key] = out[key].where(out[key] != "", vals) if key in out else vals
    frame = pd.DataFrame(out, index=df.index)
    return frame.loc[:, (frame != "").any(axis=0)]


@dataclass
class Sparc:
    samples: pd.DataFrame | None = None
    subjects: pd.DataFrame | None = None
    files: dict[str, str] = field(default_factory=dict)
    sample_columns: list[str] = field(default_factory=list)     # output names, in order
    subject_columns: dict[str, str] = field(default_factory=dict)   # subjects.xlsx column -> output name

    @property
    def columns(self) -> list[str]:
        return ["subject_id"] + [c for c in self.sample_columns if c != "subject_id"] + \
            [c for c in self.subject_columns.values() if c != "subject_id"] if self else []

    def __bool__(self) -> bool:
        return self.samples is not None or self.subjects is not None


def read_sparc(indir: InputDir, paths: dict[str, Path]) -> Sparc:
    sp = Sparc(files={k: indir.rel(p) for k, p in paths.items()})
    if "samples" in paths:
        df = _snake_frame(read_table(indir, paths["samples"]))
        if "sample_id" not in df.columns:
            raise ValueError(f"{paths['samples'].name} has no 'sample id' column")
        sp.samples = df
        sp.sample_columns = [c for c in df.columns if c != "sample_id"]
    if "subjects" in paths:
        df = _snake_frame(read_table(indir, paths["subjects"]))
        if "subject_id" not in df.columns:
            raise ValueError(f"{paths['subjects'].name} has no 'subject id' column")
        sp.subjects = df
        taken = set(sp.sample_columns)
        sp.subject_columns = {c: (f"subject_{c}" if c in taken and c not in KEYS else c) for c in df.columns}
    return sp


def _nearest(parts: tuple[str, ...], prefix: str) -> str | None:
    for p in reversed(parts):
        if p.startswith(prefix):
            return p
    return None


def join(sp: Sparc, source_rel: str, sample_id: str) -> tuple[dict[str, str], dict | None, list[str]]:
    """The SPARC columns for one sample (output names -> values), the match record for prep.json, warnings."""
    if not sp:
        return {}, None, []
    parts = Path(source_rel).parts
    sam, sub = _nearest(parts, "sam-"), _nearest(parts, "sub-")
    out: dict[str, str] = {}
    match: dict = {}
    warnings: list[str] = []
    if sp.samples is not None:
        key, how = (sam, "folder") if sam else (sample_id, "sample_id")
        rows = sp.samples[sp.samples["sample_id"] == key]
        if len(rows) > 1 and sub and "subject_id" in rows.columns:
            rows = rows[rows["subject_id"] == sub]
        if len(rows) == 1:
            row = rows.iloc[0]
            out.update({c: row[c] for c in sp.sample_columns})
            match.update({"sample": key, "sample_matched_by": how})
            sub = row.get("subject_id") or sub
        elif len(rows) > 1:
            warnings.append(f"{len(rows)} rows of the SPARC samples table match; none used")
    if sp.subjects is not None and sub:
        rows = sp.subjects[sp.subjects["subject_id"] == sub]
        if len(rows) == 1:
            row = rows.iloc[0]
            out.update({name: row[c] for c, name in sp.subject_columns.items()})
            match["subject"] = sub
        else:
            warnings.append(f"{'no row' if rows.empty else f'{len(rows)} rows'} of the SPARC subjects table match")
    if not match and not warnings:
        warnings.append("not in the SPARC metadata (no matching `sam-` / `sub-` folder or sample id)")
    if sub and "subject_id" not in out:
        out["subject_id"] = sub
    return out, (match or None), warnings
