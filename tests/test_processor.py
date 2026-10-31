"""The synthetic fixture through the processor (both input layouts), the AnnData contract, the SPARC join,
the sample sheet, the droplets, the Seurat reader and the parameter checks."""
from __future__ import annotations

import csv
import json
import logging
from collections import Counter
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from processor import handler
from processor import main as app
from processor import readers as R
from processor import inputs as I
from processor.inputs import Candidate
from processor.safe_inputs import InputDir
from processor.sparc import snake
from tests.make_fixture import DATA, HUMAN, N_PEAKS, make_fixture, write_mtx

SAMPLES = {"sam-A1", "sam-B1", "arc1", "run3", "aggr", "tab", "m1", "m2", "s1", "s2"}
VERSION = (Path(__file__).resolve().parent.parent / "VERSION").read_text().strip()


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, value)
    return app.main()


def rows(path: Path) -> dict[str, dict[str, str]]:
    with open(path, newline="") as fh:
        return {r["sample_id"]: r for r in csv.DictReader(fh)}


@pytest.fixture(scope="session")
def outputs(fixture_root, tmp_path_factory) -> dict[str, Path]:
    out = {}
    for layout in ("linked", "flat"):
        d = tmp_path_factory.mktemp(f"out-{layout}")
        mp = pytest.MonkeyPatch()
        try:
            assert run(fixture_root[layout], d, mp) == 0
        finally:
            mp.undo()
        out[layout] = d
    return out


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_every_sample_written(outputs, layout):
    out = outputs[layout]
    assert {p.name for p in (out / "cells").iterdir()} == {f"{s}.h5ad" for s in SAMPLES} | {"provenance.json"}
    assert {p.name for p in (out / "droplets").iterdir()} == \
        {f"{s}.raw.h5" for s in ("sam-A1", "sam-B1", "arc1", "run3")} | {"provenance.json"}
    assert {p.name for p in (out / "metadata").iterdir()} == \
        {"samples.csv", "samples_full.csv", "source_columns.csv", "prep.json", "provenance.json", "summary.json"}
    table = rows(out / "metadata" / "samples.csv")
    assert set(table) == SAMPLES
    prefix = "counts/" if layout == "linked" else ""
    assert {s: (r["source"], r["organism"]) for s, r in table.items()} == {
        "sam-A1": ("cellranger", "human"), "sam-B1": ("cellranger", "human"), "arc1": ("cellranger-arc", "human"),
        "run3": ("cellranger", "mouse"), "aggr": ("cellranger", "human"), "tab": ("h5ad", "human"),
        "m1": ("h5ad", "human"), "m2": ("h5ad", "human"), "s1": ("rds", "human"), "s2": ("rds", "human")}
    assert table["sam-A1"]["path"] == f"{prefix}primary/sub-A/sam-A1/count"
    assert table["arc1"]["n_genes"] == str(len(HUMAN))            # Peaks dropped

    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert prep["app"] == "sc-prep" and prep["version"] == VERSION
    assert [s["path"] for s in prep["skipped"]] == [f"{prefix}matrix_only.rds"]
    assert "not a Seurat object" in prep["skipped"][0]["reason"]
    reasons = {i["path"]: i["reason"] for i in prep["ignored"]}
    assert "raw matrix without a filtered matrix" in reasons[f"{prefix}primary/sub-B/sam-B2/count/sample_raw_feature_bc_matrix.h5"]
    assert f"{prefix}notes.txt" in reasons
    arc = prep["samples"]["arc1"]
    assert arc["atac_present"] and arc["chemistry"] == "multiome" and arc["features"]["dropped"] == {"Peaks": N_PEAKS}
    assert arc["atac_fragments"] == f"{prefix}arc1/outs/atac_fragments.tsv.gz"
    assert prep["samples"]["sam-A1"]["chemistry"] == "flex"
    assert prep["samples"]["tab"]["counts_from"] == "layers/counts"
    assert prep["samples"]["m1"]["counts_from"] == "raw/X" and prep["samples"]["m1"]["split_by"] == "sample_id"
    assert prep["samples"]["s1"]["seurat"]["layers"] == ["counts.s1", "counts.s2"]

    prov = json.loads((out / "cells" / "provenance.json").read_text())
    assert set(prov) - {"cgroup_peak_mb"} == {"app", "version", "tools", "parameters", "inputs", "started", "finished",
                                              "samples", "peak_rss_mb", "threads", "citations"}
    assert prov["app"] == "sc-prep" and prov["version"] == VERSION and prov["tools"]["scanpy"]
    assert isinstance(prov["threads"], int) and prov["peak_rss_mb"] > 0
    assert f"{prefix}primary/sub-A/sam-A1/count" in prov["inputs"] and f"{prefix}samples.xlsx" in prov["inputs"]
    assert Counter(s["status"] for s in prov["samples"]) == {"ok": len(SAMPLES), "failed": 1}
    assert prov["started"].endswith("+00:00") and len(prov["started"]) == 25
    assert {"doi": "10.1186/s13059-017-1382-0", "pmid": "29409532", "text": "Wolf FA et al. Genome Biol 2018"} \
        in prov["citations"]
    assert "10.21105/joss.07540" in [c["doi"] for c in prov["citations"]]     # a Seurat file was read
    assert all(set(c) <= {"doi", "pmid", "text"} and c["doi"] and c["text"] for c in prov["citations"])
    summary = json.loads((out / "metadata" / "summary.json").read_text())
    assert summary["n_samples"] == len(SAMPLES) and summary["samples_per_condition"]["acute"] == 1


@pytest.mark.parametrize("sid", sorted(SAMPLES))
def test_contract(outputs, sid):
    a = ad.read_h5ad(outputs["linked"] / "cells" / f"{sid}.h5ad")
    assert sp.issparse(a.X) and a.X.format == "csr" and a.X.dtype == np.int32
    assert all(o.startswith(f"{sid}:") for o in a.obs_names) and a.obs_names.is_unique
    assert list(a.obs.columns[:5]) == ["sample_id", "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
    assert isinstance(a.obs["sample_id"].dtype, pd.CategoricalDtype) and set(a.obs["sample_id"]) == {sid}
    assert np.allclose(a.obs["total_counts"], np.asarray(a.X.sum(axis=1)).ravel())
    assert list(a.var.columns) == ["gene_ids", "feature_types", "mt", "ribo"] and a.var_names.is_unique
    assert set(a.var["feature_types"]) == {"Gene Expression"}
    assert a.var["mt"].sum() >= 1 and a.var["ribo"].sum() >= 1
    pen = a.uns["pennsieve"]
    assert pen["sample_id"] == sid and pen["prep_version"] == VERSION and pen["organism"] in ("human", "mouse")


def test_sparc_join(outputs):
    out = outputs["linked"]
    a = ad.read_h5ad(out / "cells" / "sam-B1.h5ad")
    on_cells = ["condition", "subject_id", "sample_anatomical_location", "laterality", "subject_experimental_group",
                "age", "sex", "species"]
    assert list(a.obs.columns[5:]) == on_cells
    assert a.obs["condition"].iloc[0] == "chronic" and a.obs["subject_id"].iloc[0] == "sub-B"
    assert a.uns["pennsieve"]["chemistry"] == "flex"
    assert "left_average_pain" not in a.obs                      # clinical columns stay out of the cells
    table = rows(out / "metadata" / "samples.csv")
    assert list(table["sam-A1"])[6:] == on_cells                 # empty sample_experimental_group dropped
    assert table["aggr"]["condition"] == ""
    full = rows(out / "metadata" / "samples_full.csv")["sam-B1"]
    assert full["left_average_pain"] == "2" and full["freezing_method"] == "liquid nitrogen"    # coalesced
    assert full["laboratory_internal_id"] == "L2" and full["subject_laboratory_internal_id"] == "P2"
    assert "metadata_only" not in full and "sample_experimental_group" not in full
    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert prep["samples"]["sam-A1"]["sparc"] == {"sample": "sam-A1", "sample_matched_by": "folder", "subject": "sub-A"}
    assert prep["condition_from"] == "subject_experimental_group"


def test_droplets(outputs, fixture_root):
    d = outputs["linked"] / "droplets"
    src = fixture_root["linked"] / "counts"
    # Flex: copied as it is
    assert (d / "sam-A1.raw.h5").read_bytes() == (src / "primary/sub-A/sam-A1/count/sample_raw_feature_bc_matrix.h5").read_bytes()
    # ARC: Peaks left out, gene counts unchanged
    out, inp = InputDir(d), InputDir(src)
    X, bc, var, _ = R.read_10x_h5(out, d / "arc1.raw.h5")
    X0, bc0, var0, _ = R.read_10x_h5(inp, src / "arc1/outs/raw_feature_bc_matrix.h5")
    gex = (var0["feature_types"] == "Gene Expression").to_numpy()
    assert set(var["feature_types"]) == {"Gene Expression"} and list(bc) == list(bc0)
    assert (X != X0[:, np.flatnonzero(gex)]).nnz == 0
    # MTX: converted
    X, bc, var, _ = R.read_10x_h5(out, d / "run3.raw.h5")
    X0, bc0, var0, _ = R.read_10x_mtx(inp, src / "run3/outs/raw_feature_bc_matrix")
    assert (X != X0).nnz == 0 and list(bc) == list(bc0) and list(var.index) == list(var0.index)


def test_droplets_stream_in_blocks(tmp_path, monkeypatch):
    """The ARC filter gives the same matrix whatever the block size."""
    root = make_fixture(tmp_path / "fx", "flat")["flat"]
    src, indir, out = root / "arc1/outs/raw_feature_bc_matrix.h5", InputDir(root), tmp_path / "out"
    out.mkdir()
    R.write_droplets(indir, src, out / "one.h5")
    monkeypatch.setattr(R, "CHUNK_BARCODES", 7)
    R.write_droplets(indir, src, out / "blocks.h5")
    a, b = R.read_10x_h5(InputDir(out), out / "one.h5")[0], R.read_10x_h5(InputDir(out), out / "blocks.h5")[0]
    assert (a != b).nnz == 0


def test_sample_sheet_and_suffix_split(fixture_root, tmp_path, monkeypatch):
    root = tmp_path / "in"
    (root).mkdir()
    for p in (fixture_root["flat"]).iterdir():
        (root / p.name).symlink_to(p)
    pd.DataFrame({"sample_id": ["lib_a", "lib b", "drg1", "ghost"],
                  "path": ["aggr", "aggr", "primary/sub-A/sam-A1", "nowhere"],
                  "barcode_suffix": ["1", "2", "", ""],
                  "condition": ["ctrl", "pain", "pain", "ctrl"]}).to_csv(tmp_path / "samples.csv", index=False)
    (root / "samples.csv").symlink_to(tmp_path / "samples.csv")
    out = tmp_path / "out"
    assert run(root, out, monkeypatch, MAX_CELLS_PER_SAMPLE="9", MIN_COUNTS_PER_CELL="1") == 0
    table = rows(out / "metadata" / "samples.csv")
    assert set(table) == {"lib_a", "lib_b", "drg1"}               # the sheet defines the samples
    assert {s: r["condition"] for s, r in table.items()} == {"lib_a": "ctrl", "lib_b": "pain", "drg1": "pain"}
    assert table["lib_a"]["n_obs"] == "9" and table["lib_b"]["n_obs"] == "9"
    a = ad.read_h5ad(out / "cells" / "lib_a.h5ad")
    assert all(o.endswith("-1") for o in a.obs_names) and a.uns["pennsieve"]["subsampled_from"] == 12
    assert ad.read_h5ad(out / "cells" / "lib_b.h5ad").obs_names[0].endswith("-2")
    drg = ad.read_h5ad(out / "cells" / "drg1.h5ad")
    assert drg.obs["subject_id"].iloc[0] == "sub-A"              # SPARC still joined through the folders
    assert drg.obs["condition"].iloc[0] == "pain"                 # the sheet's column wins
    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert prep["samples"]["drg1"]["cells"]["n_below_min_counts"] == 1
    assert prep["samples"]["drg1"]["notes"] == []
    assert [s["sample_id"] for s in prep["skipped"]] == ["ghost"]
    assert prep["samples"]["lib_b"]["barcode_suffix"] == "2"
    assert "sample_id 'lib b' sanitised to 'lib_b'" in prep["samples"]["lib_b"]["notes"]
    assert any(i["reason"].endswith("not in the sample sheet") for i in prep["ignored"])
    # the subsample is reproducible and depends on SEED
    again = tmp_path / "again"
    assert run(root, again, monkeypatch, MAX_CELLS_PER_SAMPLE="9", MIN_COUNTS_PER_CELL="1") == 0
    assert list(ad.read_h5ad(again / "cells" / "lib_a.h5ad").obs_names) == list(a.obs_names)
    other = tmp_path / "other"
    assert run(root, other, monkeypatch, MAX_CELLS_PER_SAMPLE="9", SEED="7") == 0
    assert list(ad.read_h5ad(other / "cells" / "lib_a.h5ad").obs_names) != list(a.obs_names)


def test_seurat_v3(tmp_path):
    c = Candidate("rds", DATA / "seurat_v3.rds", {"rds": DATA / "seurat_v3.rds"}, "seurat_v3")
    reads = R.read_rds(c, InputDir(DATA), {})
    assert [r.sub_id for r in reads] == ["s1", "s2"] and reads[0].meta["seurat"]["assay_class"] == "Assay"
    assert reads[0].adata.shape == (10, 8) and "donor" in reads[0].adata.obs and "orig.ident" not in reads[0].adata.obs
    v5 = R.read_rds(Candidate("rds", DATA / "seurat_v5_split.rds", {"rds": DATA / "seurat_v5_split.rds"}, "x"),
                    InputDir(DATA), {})
    for a, b in zip(reads, v5):
        assert (a.adata.X != b.adata.X).nnz == 0 and list(a.adata.obs_names) == list(b.adata.obs_names)


def test_conditions_and_metadata_columns(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "o1"
    assert run(fixture_root["flat"], out, monkeypatch, CONDITION_COLUMN="sex", METADATA_COLUMNS="none") == 0
    a = ad.read_h5ad(out / "cells" / "sam-A1.h5ad")
    assert list(a.obs.columns[5:]) == ["condition"] and a.obs["condition"].iloc[0] == "male"
    out = tmp_path / "o2"
    assert run(fixture_root["flat"], out, monkeypatch, CONDITION_COLUMN="none",
               METADATA_COLUMNS="Left Average Pain, sex", ORGANISM="mouse") == 0
    a = ad.read_h5ad(out / "cells" / "sam-A1.h5ad")
    assert list(a.obs.columns[5:]) == ["sex", "left_average_pain"] and a.uns["pennsieve"]["organism"] == "mouse"


@pytest.mark.parametrize("params", [{"ORGANISM": "fish"}, {"MAX_CELLS_PER_SAMPLE": "-1"}, {"SEED": "x"},
                                    {"CONDITION_COLUMN": "no_such_column"}, {"SAMPLE_SHEET": "../x.csv"}])
def test_bad_parameters_exit_2(fixture_root, tmp_path, monkeypatch, params):
    assert run(fixture_root["flat"], tmp_path / "out", monkeypatch, **params) == 2


def test_no_input_exits_1(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "readme.txt").write_text("nothing here")
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_symlink_out_of_the_input_is_not_followed(tmp_path, monkeypatch, caplog):
    outside = tmp_path / "outside"
    write_mtx(outside / "filtered_feature_bc_matrix", np.ones((3, 2), dtype=np.int32), ["A-1", "B-1", "C-1"],
              ["G1", "G2"], ["ENSG1", "ENSG2"], ["Gene Expression"] * 2)
    inp = tmp_path / "in" / "data"
    inp.mkdir(parents=True)
    write_mtx(inp / "s1" / "filtered_feature_bc_matrix", np.ones((2, 2), dtype=np.int32), ["A-1", "B-1"],
              ["G1", "G2"], ["ENSG1", "ENSG2"], ["Gene Expression"] * 2)
    (inp / "deep").mkdir()
    os.symlink(outside, inp / "deep" / "sneaky")          # below the port folder: must stay inside the input
    out = tmp_path / "out"
    with caplog.at_level(logging.WARNING):
        assert run(tmp_path / "in", out, monkeypatch) == 0
    assert set(rows(out / "metadata" / "samples.csv")) == {"s1"}
    assert "skipped data/deep/sneaky: resolves outside the input" in caplog.text


def test_non_count_h5ad_fails_the_sample(tmp_path, monkeypatch):
    inp = tmp_path / "in"
    inp.mkdir()
    a = ad.AnnData(X=sp.csr_matrix(np.full((4, 3), 0.5, dtype=np.float32)), var=pd.DataFrame(index=["A", "B", "C"]))
    a.write_h5ad(inp / "norm.h5ad")
    assert run(inp, tmp_path / "out", monkeypatch) == 1


def test_snake():
    assert snake("subject experimental group") == "subject_experimental_group"
    assert snake("wasDerivedFrom") == "was_derived_from"
    assert snake("RRID for strain") == "rrid_for_strain"
    assert snake(" Freezing method ") == "freezing_method"


def test_lambda_handler(fixture_root, tmp_path, monkeypatch):
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("METADATA_COLUMNS", "sex")      # left from a previous invocation: must be cleared
    res = handler.handler({"inputDir": str(fixture_root["flat"]), "outputDir": str(tmp_path / "out"),
                           "MAX_CELLS_PER_SAMPLE": 5}, None)
    assert res["status"] == "ok"
    a = ad.read_h5ad(tmp_path / "out" / "cells" / "sam-A1.h5ad")
    assert a.n_obs == 5 and "age" in a.obs
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(tmp_path / "missing"), "outputDir": str(tmp_path / "o")}, None)


def test_h5_written_is_10x_format(outputs):
    with h5py.File(outputs["linked"] / "droplets" / "run3.raw.h5") as h:
        assert set(h["matrix"]) == {"barcodes", "data", "indices", "indptr", "shape", "features"}
        assert set(h["matrix/features"]) == {"name", "id", "feature_type", "genome", "_all_tag_keys"}


def test_source_columns(outputs, fixture_root, tmp_path, monkeypatch):
    """Label-like source columns are copied; others (clinical fields) only when named in SOURCE_COLUMNS, and
    their names (never values) are listed in source_columns.csv."""
    out = outputs["linked"]
    tab = ad.read_h5ad(out / "cells" / "tab.h5ad")
    assert "source_cell_type" in tab.obs and "pain_score" not in tab.obs
    s1 = ad.read_h5ad(out / "cells" / "s1.h5ad")
    assert "donor" not in s1.obs and "nCount_RNA" not in s1.obs
    listed = pd.read_csv(out / "metadata" / "source_columns.csv", keep_default_na=False)
    assert list(listed.columns) == ["sample_id", "column", "copied_as"]
    by = {(r.sample_id, r.column): r.copied_as for r in listed.itertuples()}
    assert by[("tab", "cell_type")] == "source_cell_type" and by[("tab", "pain_score")] == ""
    assert by[("s1", "donor")] == ""
    assert "7" not in (out / "metadata" / "source_columns.csv").read_text()   # no values
    o2 = tmp_path / "o2"
    assert run(fixture_root["flat"], o2, monkeypatch, SOURCE_COLUMNS="pain_score, donor") == 0
    assert ad.read_h5ad(o2 / "cells" / "tab.h5ad").obs["pain_score"].iloc[0] == 7.0
    assert "donor" in ad.read_h5ad(o2 / "cells" / "s1.h5ad").obs


@pytest.mark.parametrize("sheet", ["oversize", "binary", "no_path"])
def test_unusable_sample_sheet_exits_cleanly(fixture_root, tmp_path, monkeypatch, sheet):
    root = tmp_path / "in"
    root.mkdir()
    (root / "sam").symlink_to(fixture_root["flat"] / "primary/sub-A/sam-A1")
    if sheet == "oversize":
        monkeypatch.setattr(I, "MAX_TABLE_BYTES", 10)
        (root / "samples.csv").write_text("sample_id,path\nx,sam\n")
    elif sheet == "binary":
        (root / "samples.csv").write_bytes(b"\xff\xfe\x00\x81bad")
    else:                                   # a sheet that is not named samples.csv cannot be SPARC metadata
        (root / "sheet.csv").write_text("sample_id,where\nx,sam\n")
    params = {"SAMPLE_SHEET": "sheet.csv"} if sheet == "no_path" else {}
    assert run(root, tmp_path / "out", monkeypatch, **params) == 1


def test_h5ad_with_external_link_is_refused(tmp_path, monkeypatch):
    inp = tmp_path / "in"
    inp.mkdir()
    (tmp_path / "secret.h5").write_bytes(b"")
    a = ad.AnnData(X=sp.csr_matrix(np.ones((3, 2), dtype=np.int32)), var=pd.DataFrame(index=["A", "B"]))
    a.write_h5ad(inp / "ok.h5ad")
    a.write_h5ad(inp / "linked.h5ad")
    with h5py.File(inp / "linked.h5ad", "a") as h:
        h["uns/escape"] = h5py.ExternalLink(str(tmp_path / "secret.h5"), "/")
    out = tmp_path / "out"
    assert run(inp, out, monkeypatch) == 0
    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert set(prep["samples"]) == {"ok"} and "external link" in prep["skipped"][0]["reason"]
