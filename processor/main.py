"""
sc-prep: the first step of the single-cell workflows.

Finds every sample under INPUT_DIR (Cell Ranger count / multi (Flex) / ARC matrices, .h5ad, Seurat
.rds), joins the optional sample sheet and the SPARC samples / subjects sheets, and writes

  OUTPUT_DIR/cells/<sample_id>.h5ad        raw counts (CSR int32), QC metrics, metadata columns, uns["pennsieve"]
  OUTPUT_DIR/droplets/<sample_id>.raw.h5   the raw (unfiltered) matrix as 10x HDF5, when the sample has one
  OUTPUT_DIR/metadata/samples.csv          one row per sample: required columns + the metadata put on cells
  OUTPUT_DIR/metadata/samples_full.csv     the same with every sample sheet and SPARC column
  OUTPUT_DIR/metadata/prep.json            per sample: what was read, kept, dropped; skipped and ignored inputs
  OUTPUT_DIR/<port>/provenance.json, summary.json

One sample at a time; nothing is normalised, and cells are only removed by MIN_COUNTS_PER_CELL
and MAX_CELLS_PER_SAMPLE.
"""
from __future__ import annotations

from processor import threads  # noqa: F401  (sets the thread environment before numpy is imported)

import csv  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import os  # noqa: E402
import platform as _platform  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from collections import Counter  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from importlib import metadata as importlib_metadata  # noqa: E402
from pathlib import Path  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-sc-prep")
# readers build AnnData with the source's (possibly repeated) symbols; contract.finalise makes them unique
warnings.filterwarnings("ignore", message="Variable names are not unique")
for noisy in ("numba", "matplotlib", "h5py", "anndata", "scanpy", "openpyxl"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "sc-prep"
REPO_ROOT = Path(__file__).resolve().parent.parent
ORGANISMS = ("auto", "human", "mouse", "rat")
CONDITION_AUTO = ("condition", "subject_experimental_group", "sample_experimental_group")
CITATIONS = [
    {"doi": "10.1186/s13059-017-1382-0", "pmid": "29409532", "text": "Wolf FA et al. Genome Biol 2018"},
    {"doi": "10.21105/joss.04371", "text": "Virshup I et al. J Open Source Softw 2024"},
    {"doi": "10.1038/s41576-023-00586-w", "pmid": "37002403", "text": "Heumos L et al. Nat Rev Genet 2023"},
]
RDS_CITATION = {"doi": "10.21105/joss.07540", "text": "Ramos-Carreño C et al. J Open Source Softw 2024"}

DEFAULTS = {
    "SAMPLE_SHEET": "samples.csv",
    "SAMPLE_SUFFIX_COLUMN": "barcode_suffix",
    "METADATA_COLUMNS": "subject_id,subject_experimental_group,sample_experimental_group,sex,age,species,"
                        "sample_anatomical_location,laterality",
    "CONDITION_COLUMN": "auto",
    "SOURCE_COLUMNS": "",
    "MAX_CELLS_PER_SAMPLE": "0",
    "MIN_COUNTS_PER_CELL": "0",
    "ORGANISM": "auto",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    from processor.sparc import snake

    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p_source = [c.strip() for c in raw["SOURCE_COLUMNS"].split(",") if c.strip()]
    p: dict = {k: raw[k].strip() for k in ("SAMPLE_SHEET", "SAMPLE_SUFFIX_COLUMN", "CONDITION_COLUMN")}
    if "/" in p["SAMPLE_SHEET"] or p["SAMPLE_SHEET"] in (".", ".."):
        raise ParameterError(f"SAMPLE_SHEET is a file name, not a path: {p['SAMPLE_SHEET']!r}")
    cols = raw["METADATA_COLUMNS"].strip()
    p["METADATA_COLUMNS"] = [] if cols.lower() == "none" else list(dict.fromkeys(snake(c) for c in cols.split(",") if c.strip()))
    p["SOURCE_COLUMNS"] = list(dict.fromkeys(p_source))
    p["ORGANISM"] = raw["ORGANISM"].strip().lower()
    if p["ORGANISM"] not in ORGANISMS:
        raise ParameterError(f"ORGANISM must be one of {', '.join(ORGANISMS)}, got {raw['ORGANISM']!r}")
    try:
        for k in ("MAX_CELLS_PER_SAMPLE", "MIN_COUNTS_PER_CELL", "SEED", "THREADS"):
            p[k] = int(raw[k])
    except ValueError as exc:
        raise ParameterError(f"bad integer parameter: {exc}") from exc
    if min(p["MAX_CELLS_PER_SAMPLE"], p["MIN_COUNTS_PER_CELL"], p["THREADS"]) < 0:
        raise ParameterError("MAX_CELLS_PER_SAMPLE, MIN_COUNTS_PER_CELL and THREADS must be >= 0")
    return p


def version_of(root: Path) -> str:
    try:
        return (root / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict[str, str]:
    out = {"python": _platform.python_version()}
    for pkg in ("scanpy", "anndata", "numpy", "scipy", "pandas", "h5py", "openpyxl", "defusedxml", "rdata"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def peak_rss_mb() -> float:
    """The larger peak RSS of this process and its children (the Seurat reader)."""
    import resource

    peak = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)


def cgroup_peak_mb() -> float | None:
    """The container's peak memory (cgroup v2 memory.peak, v1 max_usage_in_bytes; includes page cache)."""
    for f in ("/sys/fs/cgroup/memory.peak", "/sys/fs/cgroup/memory/memory.max_usage_in_bytes"):
        try:
            return round(int(open(f).read()) / 2**20, 1)
        except (OSError, ValueError):
            pass
    return None


def write_csv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
        fh.write("\n")


def resolve_condition(param: str, available: list[str]) -> str | None:
    """The column copied to `condition`: CONDITION_COLUMN, or for auto the first of CONDITION_AUTO present."""
    from processor.sparc import snake

    if param.lower() == "none":
        return None
    if param.lower() == "auto":
        return next((c for c in CONDITION_AUTO if c in available), None)
    for name in (param, snake(param)):
        if name in available:
            return name
    raise ParameterError(f"CONDITION_COLUMN {param!r} is not a column of the sample sheet or SPARC metadata "
                         f"(columns: {', '.join(available) or 'none'})")


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir or not Path(input_dir).is_dir():
        log.error("INPUT_DIR (an existing folder) and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    n_threads = threads.THREADS     # fixed at start-up, before numpy was imported (see processor/threads.py)

    from processor import contract as C
    from processor import inputs as I
    from processor import readers as R
    from processor import sparc as S
    from processor.safe_inputs import InputDir

    app_version = version_of(REPO_ROOT)
    indir = InputDir(input_dir)
    rel = indir.rel

    # 1. what is there
    found = I.walk(indir, params["SAMPLE_SHEET"])
    log.info("input: %d files, %d input(s) detected (%s), sample sheet: %s, SPARC: %s, threads %d",
             found.n_files, len(found.candidates),
             ", ".join(f"{k} x{v}" for k, v in sorted(Counter(c.source for c in found.candidates).items())) or "none",
             rel(found.sheet) if found.sheet else "none",
             ", ".join(rel(p) for p in found.sparc.values()) or "none", n_threads)
    if not found.candidates:
        log.error("no single-cell sample found under INPUT_DIR (%d files; e.g. %s)", found.n_files,
                  "; ".join(f"{i['path']}: {i['reason']}" for i in found.ignored[:10]) or "nothing")
        return 1
    try:
        # a SPARC samples.csv (no `path` column) is SPARC metadata, not the sample sheet
        if found.sheet is not None and found.sheet == found.sparc.get("samples"):
            if "path" in I.read_table(indir, found.sheet).columns:
                del found.sparc["samples"]
            else:
                found.sheet = None
        sheet = I.read_sheet(indir, found.sheet, params["SAMPLE_SUFFIX_COLUMN"])
        sparc = S.read_sparc(indir, found.sparc)
    except (ValueError, OSError, RuntimeError) as exc:     # UnsafeInput, UnicodeDecodeError are ValueErrors
        log.error("cannot use a metadata table: %s", exc)
        return 1
    sheet_errors = I.assign_identities(found.candidates, sheet, indir, params["SAMPLE_SUFFIX_COLUMN"])
    skipped: list[dict] = []
    for e in sheet_errors:
        log.error("sample sheet row %s (path %s): %s", e["sample_id"] or "?", e["path"], e["reason"])
        skipped.append(dict(e))
    if sheet.path:
        for c in found.candidates:
            if c.sheet_row is None and not c.suffix_rows:
                log.warning("%s: not in %s; ignored", rel(c.path), sheet.path.name)
                found.ignored.append({"path": rel(c.path), "reason": f"{c.source} input not in the sample sheet"})
                found.n_ignored += 1
        found.candidates = [c for c in found.candidates if c.sheet_row is not None or c.suffix_rows]

    # 2. which metadata columns go where
    available = list(dict.fromkeys(sheet.columns + sparc.columns))
    try:
        condition_from = resolve_condition(params["CONDITION_COLUMN"], available)
    except ParameterError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    on_cells = sheet.columns + [c for c in sparc.columns if c in params["METADATA_COLUMNS"] and c not in sheet.columns]
    obs_columns = list(dict.fromkeys((["condition"] if condition_from else []) + on_cells))
    log.info("metadata: %d sample sheet column(s), %d SPARC column(s), %d put on cells, condition from %s",
             len(sheet.columns), len(sparc.columns), len(obs_columns), condition_from or "nothing")

    def metadata_for(sid: str, source_rel: str, row: dict | None) -> tuple[dict, dict | None, list[str]]:
        values, match, warnings = S.join(sparc, source_rel, sid)
        values.update({k: v for k, v in (row or {}).items() if k in sheet.columns and v})
        if condition_from:
            values["condition"] = values.get(condition_from, "")
        return values, match, warnings

    cells_dir, drop_dir, meta_dir = (output_dir / n for n in ("cells", "droplets", "metadata"))
    source_rows: list[dict] = []
    for d in (cells_dir, drop_dir, meta_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 3. one sample at a time
    rows: list[dict] = []
    prep_samples: dict[str, dict] = {}
    prov_samples: list[dict] = []
    written: set[str] = set()

    def fail(sid: str, src: str, source: str, exc: Exception, t0: float) -> None:
        log.error("%s (%s): FAILED: %s", sid, src, exc)
        skipped.append({"sample_id": sid, "path": src, "reason": f"{type(exc).__name__}: {exc}"})
        prov_samples.append({"sample_id": sid, "input": src, "source": source, "status": "failed",
                             "error": f"{type(exc).__name__}: {exc}", "seconds": round(time.monotonic() - t0, 2)})

    for c in found.candidates:
        t0 = time.monotonic()
        src = rel(c.path)
        for note in c.notes:
            log.warning("%s: %s", src, note)
        try:
            reads = R.READERS[c.source](c, indir, params)
        except Exception as exc:  # noqa: BLE001 - one sample's failure must not sink the run
            fail(c.sample_id, src, c.source, exc, t0)
            continue
        for r in reads:
            t1 = t0 if len(reads) == 1 else time.monotonic()   # a split input's read time is not attributed
            sid = c.sample_id if r.sub_id is None else r.sub_id
            row = c.suffix_rows.get(r.meta.get("barcode_suffix")) if r.meta.get("barcode_suffix") else c.sheet_row
            if sid in written:
                fail(sid, src, c.source, ValueError(f"sample_id {sid} already written"), t1)
                continue
            try:
                values, match, sparc_warnings = metadata_for(sid, src if r.sub_id is None else "", row)
                adata, info = C.finalise(r, sid, src, params, values, obs_columns, values.get("species", ""), app_version)
                adata.write_h5ad(cells_dir / f"{sid}.h5ad", compression="gzip")
                drop = None
                if "raw" in c.files and r.sub_id is None:
                    drop = R.write_droplets(indir, c.files["raw"], drop_dir / f"{sid}.raw.h5")
                    drop.update({"from": rel(c.files["raw"]), "output": f"droplets/{sid}.raw.h5"})
            except Exception as exc:  # noqa: BLE001
                fail(sid, src, c.source, exc, t1)
                continue
            written.add(sid)
            info["warnings"] += sparc_warnings
            if drop is None:
                info["warnings"].append("no raw matrix: sc-qc skips ambient RNA removal for this sample")
            seconds = round(time.monotonic() - t1, 2)
            source_rows += [{"sample_id": sid, "column": col, "copied_as": name or ""}
                            for col, name in info["source_columns"].items()]
            info.update({"input": src, "output": f"cells/{sid}.h5ad", "droplets": drop, "sparc": match,
                         "in_sample_sheet": (row is not None) if sheet.path else None, "notes": list(c.notes),
                         "seconds": seconds})
            prep_samples[sid] = info
            for w in info["warnings"]:
                log.warning("%s: %s", sid, w)
            rows.append({"sample_id": sid, "path": src, "source": info["source"], "n_obs": info["n_obs"],
                         "n_genes": info["n_genes"], "organism": info["organism"], **values})
            prov_samples.append({"sample_id": sid, "input": src, "source": info["source"], "status": "ok",
                                 "n_obs": info["n_obs"], "n_genes": info["n_genes"], "droplets": drop is not None,
                                 "seconds": seconds})
            cells = info["cells"]
            log.info("%s: %s %s -> %d cells x %d genes (%d source cells, %d below MIN_COUNTS_PER_CELL%s; %d non-gene "
                     "features dropped), median %d counts / %d genes / %.1f%% mito, organism %s, chemistry %s, "
                     "droplets %s, SPARC %s, %.1fs",
                     sid, info["source"], src, info["n_obs"], info["n_genes"], cells["n_source"],
                     cells["n_below_min_counts"],
                     f", subsampled to {cells['n_subsampled_to']}" if cells["n_subsampled_to"] else "",
                     sum(info["features"]["dropped"].values()), adata.obs["total_counts"].median(),
                     adata.obs["n_genes_by_counts"].median(), adata.obs["pct_counts_mt"].median(), info["organism"],
                     info["chemistry"] or "unknown", f"{drop['n_barcodes']} barcodes" if drop else "none",
                     "matched" if match else "no", seconds)
            del adata

    # 4. metadata tables
    required = list(I.REQUIRED_COLUMNS)
    filled = lambda cols: [c for c in cols if any(r.get(c) for r in rows)]  # noqa: E731
    write_csv(meta_dir / "samples.csv", required + filled(obs_columns), rows)
    full_cols = list(dict.fromkeys((["condition"] if condition_from else []) + sheet.columns + sparc.columns))
    write_csv(meta_dir / "samples_full.csv", required + filled(full_cols), rows)
    write_csv(meta_dir / "source_columns.csv", ["sample_id", "column", "copied_as"], source_rows)
    finished = datetime.now(timezone.utc)
    common = {"app": APP_NAME, "version": app_version, "finished": finished.isoformat(timespec="seconds")}
    write_json(meta_dir / "prep.json", {
        **common, "parameters": params, "sample_sheet": rel(found.sheet) if found.sheet else None, "sample_sheet_columns": sheet.columns,
        "sparc": sparc.files, "sparc_columns": sparc.columns, "condition_from": condition_from,
        "columns_on_cells": obs_columns, "samples": prep_samples, "skipped": skipped,
        "ignored": found.ignored, "n_ignored": found.n_ignored})
    conditions = Counter(r.get("condition", "") for r in rows) if condition_from else Counter()
    summary = {**common, "n_samples": len(rows), "n_failed": len(skipped),
               "n_cells": int(sum(r["n_obs"] for r in rows)),
               "n_with_droplets": sum(1 for s in prep_samples.values() if s["droplets"]),
               "sources": dict(Counter(r["source"] for r in rows)),
               "organisms": dict(Counter(r["organism"] for r in rows)),
               "condition_from": condition_from, "samples_per_condition": dict(conditions),
               "samples": [{"sample_id": r["sample_id"], "n_obs": r["n_obs"], "n_genes": r["n_genes"],
                            "n_obs_source": prep_samples[r["sample_id"]]["cells"]["n_source"]} for r in rows]}
    # sheet rows that matched nothing are listed with the samples, status "skipped"
    prov_samples += [{"sample_id": e["sample_id"], "input": e["path"], "status": "skipped", "error": e["reason"]}
                     for e in sheet_errors]
    provenance = {
        "app": APP_NAME, "version": app_version, "tools": tool_versions(), "parameters": params,
        "inputs": [rel(c.path) for c in found.candidates] + ([rel(found.sheet)] if found.sheet else [])
        + [rel(p) for p in found.sparc.values()],
        "started": started.isoformat(timespec="seconds"), "finished": finished.isoformat(timespec="seconds"),
        "samples": prov_samples, "peak_rss_mb": peak_rss_mb(), "threads": n_threads,
        "citations": CITATIONS + ([RDS_CITATION] if any(r["source"] == "rds" for r in rows) else []),
    }
    cgroup_peak = cgroup_peak_mb()
    if cgroup_peak is not None:
        provenance["cgroup_peak_mb"] = cgroup_peak
    for d in (cells_dir, drop_dir, meta_dir):
        write_json(d / "provenance.json", provenance)
    write_json(meta_dir / "summary.json", summary)
    log.info("%d sample(s) written (%d cells, %d with droplets), %d skipped, peak RSS %s MB, %.1fs",
             len(rows), summary["n_cells"], summary["n_with_droplets"], len(skipped), provenance["peak_rss_mb"],
             time.monotonic() - t_start)
    if not rows:
        log.error("every sample failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
