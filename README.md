# sc-prep

The first step of the single-cell workflows (A. Single-nucleus Atlas Mapping & Condition
Comparison: `sc-prep → sc-qc → sc-integrate → sc-annotate → {sc-markers, st-compare} → sc-report`;
B. Multiome Gene Regulation takes its RNA half from the same chain). It takes single-cell or
single-nucleus RNA-seq data as Cell Ranger, scanpy or Seurat leave it, together with the
dataset's sample metadata, and turns every sample into one AnnData file of the same shape, so
every later step can read 10x 3'/5', Flex, Multiome, `.h5ad` and Seurat data without knowing
where it came from.

**For a biologist:** point it at a dataset folder, for example a SPARC dataset from Pennsieve
Discover with its `samples.xlsx` and `subjects.xlsx`. It finds each sample, keeps the raw
counts and adds a few quality numbers per cell (total counts, genes detected, % mitochondrial
and % ribosomal reads). It also attaches what is known about the sample and the donor: subject,
sex, age, laterality, and a `condition` column (by default the SPARC *experimental group*, e.g.
acute vs chronic pain) that the comparison steps use. Clinical fields (pain scores, medication,
comorbidities) are kept in `samples_full.csv` but **not** copied onto the cells unless you name
them in `METADATA_COLUMNS`. It also copies each sample's raw droplet matrix, which `sc-qc` needs
to remove ambient RNA, a large source of contamination in nuclei from neurons. It does not
filter or normalise anything; that is `sc-qc`'s job.

**For a computational analyst:** one `<sample_id>.h5ad` per sample following the AnnData
contract of the workflows (`docs/conventions.md` in `workflow-single-cell`): raw integer
counts as CSR int32 in `X`; only `Gene Expression` features (ARC `Peaks` are left for
`atac-prep`, other feature types are dropped and counted); unique gene symbols as `var_names`
with `gene_ids`, `mt`, `ribo`; `obs_names` `<sample_id>:<barcode>`; QC metrics from
`scanpy.pp.calculate_qc_metrics`; metadata columns as categoricals; `uns["pennsieve"]` with
`chemistry` and `organism`. Samples are read and written one at a time and matrices are never
densified. The raw matrix goes to `droplets/<sample_id>.raw.h5` (10x HDF5; copied as it is, or
with ARC's peaks streamed out in barcode blocks, or converted from MTX).

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `counts` | input | root | The data folder tree (searched recursively; symlinks followed inside the input), with the optional sample sheet and SPARC `samples` / `subjects` sheets anywhere in it. |
| `cells` | output | `OUTPUT_DIR/cells/` | `<sample_id>.h5ad` per sample; `provenance.json`, `summary.json`. Link it to `sc-qc` (and, after the later steps, to st-compare's `spatial` input). |
| `droplets` | output | `OUTPUT_DIR/droplets/` | `<sample_id>.raw.h5` per sample with a raw matrix; `provenance.json`, `summary.json`. |
| `metadata` | output | `OUTPUT_DIR/metadata/` | `samples.csv`, `samples_full.csv`, `prep.json`, `provenance.json`, `summary.json`. |

## Inputs it recognises

| Source | Detected by | Read |
|---|---|---|
| Cell Ranger count (3', 5') | a folder with `filtered_feature_bc_matrix.h5`, or the MTX folder `filtered_feature_bc_matrix/` (`matrix.mtx.gz`, `features.tsv.gz` / `genes.tsv`, `barcodes.tsv.gz`) | The filtered matrix (h5 preferred over MTX). The raw matrix next to it (`raw_feature_bc_matrix.h5` or `raw_feature_bc_matrix/`) becomes the sample's droplets. |
| Cell Ranger multi / Flex | `count/sample_filtered_feature_bc_matrix.h5` (or `count/sample_filtered_feature_bc_matrix/`) in `per_sample_outs/<sample>/` or a dataset's `sam-…/count/` | As above; raw from `sample_raw_feature_bc_matrix.h5`. Chemistry `flex` from the h5's `chemistry_description`. |
| Cell Ranger ARC (Multiome) | `filtered_feature_bc_matrix.h5` with `Gene Expression` and `Peaks` features | Gene Expression only; `atac_present` and the path of `atac_fragments.tsv.gz` recorded in `prep.json` for `atac-prep`. The raw h5 is written without its peaks. Chemistry `multiome`. |
| Aggregated matrix | any of the above whose barcodes end in `-1`, `-2`, … (`cellranger aggr`, or a dataset's combined matrix such as Discover 465's) | One sample, unless the sample sheet has rows with the same `path` and a `SAMPLE_SUFFIX_COLUMN` value: then each suffix is its own sample (barcodes with an unlisted suffix are dropped and counted). Its raw matrix is not split (no droplets). |
| AnnData | `*.h5ad` | Only `obs`, `var` and the counts are read, from the first integer-valued of `layers["counts"]`, `X`, `raw.X` (a file without integer counts fails). Ensembl `var_names` are swapped for a symbol column (`feature_name`, `gene_symbols`, …). A file with several `obs["sample_id"]` values is split. The file's own per-cell columns are kept (a column whose name a step writes, such as `cell_type`, becomes `source_cell_type`). |
| Seurat | `*.rds` holding a Seurat object (v3 / v4 `Assay`, v5 `Assay5` with one `counts` layer or split `counts.<x>` layers) | The `counts` of the `RNA` assay (else the active assay) and `meta.data`; split by `sample_id`, else `orig.ident` when they hold more than one value. See *Seurat files* below. |

Detection is by file and folder names; one sample per Cell Ranger folder (not descended into
further) or per `.h5ad` / `.rds` file. A raw matrix without a filtered one is not used (sc-prep
does not call cells). Everything not used is listed in `prep.json → ignored` with the reason.
Without a sample sheet, `sample_id` = the folder name above `outs/` / `count/` (sanitised to
`[A-Za-z0-9_-]`), or the file stem.

## Sample metadata

In order of precedence:

1. **A sample sheet** (`SAMPLE_SHEET`, default `samples.csv`; CSV, TSV or `.xlsx`) with
   `sample_id` and `path` (matches when it equals or ends the sample's path, with or without
   `outs/` / `count/`). Its other columns go onto every cell. **With a sheet, only its rows are
   samples**; anything else detected is listed in `prep.json → ignored`, and a row that matches
   nothing (or several samples) is listed in `prep.json → skipped`.
2. **SPARC dataset metadata**: `samples.xlsx` / `subjects.xlsx` (or `.csv` / `.tsv`; the
   shallowest of each name), the SPARC Dataset Structure. A sample whose path has a `sam-…`
   folder matches the `sample id` row of the nearest one (Discover 480 nests derived samples:
   `sam-UTD-DN0038/sam-UTD-DN0038-L5/count/` matches `sam-UTD-DN0038-L5`); otherwise the
   `sample_id` itself is looked up. Its `subject id`, or a `sub-…` folder above it, joins the
   subject row (if a sample id occurs twice, the `sub-…` folder decides). Column names become
   snake_case (`subject experimental group` → `subject_experimental_group`, `wasDerivedFrom` →
   `was_derived_from`); columns that fall together (`Freezing method` / `freezing method`) are
   merged, first non-empty value wins; subject columns that clash with a sample column get a
   `subject_` prefix; empty columns are dropped. **Only the `METADATA_COLUMNS` allow-list goes
   onto cells and into `samples.csv`**; every column goes to `samples_full.csv`. How each
   sample matched is in `prep.json → samples.<id>.sparc`.
3. **`condition`**: `CONDITION_COLUMN` (`auto`: the sheet's `condition`, else
   `subject_experimental_group`, else `sample_experimental_group`) is copied to a `condition`
   column, which st-compare compares by default.

Values are never logged (subject sheets hold clinical fields); logs give counts only.

## Parameters

| Name | Default | What it does |
|---|---|---|
| `SAMPLE_SHEET` | `samples.csv` | File name of the sample sheet (see above). A SPARC `samples.csv` (no `path` column) is read as SPARC metadata instead. |
| `SAMPLE_SUFFIX_COLUMN` | `barcode_suffix` | Sample sheet column that splits an aggregated matrix by barcode suffix. |
| `METADATA_COLUMNS` | `subject_id,subject_experimental_group,sample_experimental_group,sex,age,species,sample_anatomical_location,laterality` | SPARC columns put on cells and in `samples.csv` (snake_case; names given in any case or spacing are converted). `none` for none. Columns missing from the metadata are skipped; columns empty for every sample are left out. |
| `CONDITION_COLUMN` | `auto` | Column copied to `condition`; `none` for no condition column. A name that is not a column exits with code 2 and lists the columns. |
| `MAX_CELLS_PER_SAMPLE` | `0` | At most this many cells per sample: a seeded random subset, in source order (`uns["pennsieve"]["subsampled_from"]` records the count before). `0` keeps all. Droplets are never subsampled. |
| `MIN_COUNTS_PER_CELL` | `0` | Drops cells with fewer total counts, before the subsample. Quality filtering proper is sc-qc's. |
| `ORGANISM` | `auto` | `human`, `mouse` or `rat`; `auto` takes the metadata's `species`, else the 10x genome (`GRCh38`, `GRCm39`, `mRatBN7.2`, …), else the Ensembl id prefix, else upper-case symbols → human; otherwise `unknown` with a warning. |
| `SEED` | `42` | Seed of the subsample, combined with the `sample_id`. |
| `THREADS` | `0` | Threads for BLAS / numba; `0` = from the container's CPU quota (2 on Lambda). |

Bad parameters exit with code 2; no usable input, or every sample failing, exits with 1. A
sample that fails (an unreadable file, an `.h5ad` without integer counts, an `.rds` that is not
a Seurat object) is logged, listed in `prep.json → skipped`, and the run goes on.

## Outputs

**`cells/<sample_id>.h5ad`**

| Slot | Content |
|---|---|
| `X` | Raw counts, CSR int32. |
| `obs_names` | `<sample_id>:<barcode>` (an `.h5ad` this app wrote before keeps its ids). |
| `obs` | `sample_id`, `total_counts`, `n_genes_by_counts`, `pct_counts_mt`, `pct_counts_ribo`, then `condition`, the sample sheet's columns and the `METADATA_COLUMNS` (categorical; missing values NaN), then the source's own per-cell columns (`.h5ad` / Seurat). Cells in source order. |
| `var` | Gene symbols (made unique with `-1`, … suffixes; a warning counts them), `gene_ids`, `feature_types` (`Gene Expression`), `mt` (`MT-`/`mt-`), `ribo` (`RPS`/`RPL`), both case-insensitive. |
| `uns["pennsieve"]` | `source` (input path), `source_type` (`cellranger`, `cellranger-arc`, `h5ad`, `rds`), `chemistry` when known (`3p`, `5p`, `flex`, `multiome`), `organism`, `prep_version`, `sample_id`, `n_obs_source`, `subsampled_from`. |

**`droplets/<sample_id>.raw.h5`**: 10x HDF5 (format 2), Gene Expression features, every barcode
of the raw matrix.

**`metadata/samples.csv`**: `sample_id, path, source, n_obs, n_genes, organism`, then
`condition`, the sheet's columns and the allow-listed SPARC columns that are filled for at
least one sample. **`metadata/samples_full.csv`**: the same rows with every sheet and SPARC
column.

**`metadata/prep.json`**: parameters; the sheet and SPARC files used and their columns; the
condition source; per sample: source, files read, chemistry, Cell Ranger version, organism and
how it was decided, `features` (`n_source`, `n_genes`, `dropped` by type), `cells`
(`n_source`, `n_below_min_counts`, `n_subsampled_to`, `n_written`), `atac_present`,
`atac_fragments`, `droplets` (copied or converted, barcodes, features dropped), `sparc`
(sample / subject matched and how), warnings, seconds; `skipped`; `ignored` (at most 500
listed, `n_ignored` counts all).

**`provenance.json`** (each port): `app`, `version`, `tools`, `parameters`, `inputs`,
`started`, `finished`, `samples` (per-sample status and seconds), `skipped`, `peak_rss_mb`
(the largest of this process and the Seurat reader), `citations` (`{doi, pmid, text}`).
**`summary.json`** (each port): numbers of samples, cells, samples with droplets, sources,
organisms, samples per condition, and per sample `n_obs`, `n_genes`, `n_obs_source`.

## Seurat files

`.rds` files are untrusted input, and R's own `readRDS` can run code stored in a file (R
serialization can carry promises that are evaluated on access; CVE-2024-27322). sc-prep does not
use R: it parses the file with [rdata](https://github.com/vnmabus/rdata) (pure Python, MIT, a
JOSS-reviewed reader of R's serialization format) and converts only the slots it needs (the
assay's `counts` and `meta.data`; closures, environments and commands are never touched). The
conversion runs in a child process (`processor/rds.py`, run with `python -I`, an empty
environment, a memory cap of 80% of the container's memory via `RLIMIT_AS`, and a 30-minute
timeout), which writes a temporary `.h5ad` under `/tmp` that the main process reads like any
other. Tested on objects written by SeuratObject 5.4.0 (R 4.5.1): v5 single layer, v5 split
layers, v3 assay, a 20 000 × 20 000 object with reductions, closures and environments in `misc`
(parsed in 1.8 s, 1.6 GB peak), and an xz-compressed file. `SingleCellExperiment` and other R
classes are refused. A Seurat `SCT` assay is not raw counts: when an object has no `RNA` assay
sc-prep reads the active one and warns.

## Runtime

`runtime: cpu 2048, memory 3008 MiB, computeTypes [standard, lambda], timeout 3600 s`. Measured
with the linux/amd64 image under amd64 emulation (Docker Desktop on an M-series Mac, shared with
other jobs), `--memory 3008m --cpus 2 --network none --read-only`, full data (no subsampling);
time is the process's own wall clock (the host's was several times longer while Docker was
contended), peak RSS from the process, cgroup peak from `memory.peak` (includes page cache):

| Demo | Samples | Nuclei | Largest sample | Time | Peak RSS | cgroup peak |
|---|---|---|---|---|---|---|
| Discover 475 (DRG C2, 10x Flex) | 12 | 67 902 | 21 553 nuclei, 45 M non-zeros (MTX, 36 s) | 63 s | 1.65 GB | 1.77 GB |
| Discover 480 (DRG atlas, 10x Flex), the 8 samples downloaded | 8 | 72 178 | 27 745 nuclei, 78 M non-zeros (12 s; 1.85 GB peak RSS run alone) | 39 s | 1.93 GB | 2.34 GB |

Memory follows the largest single sample: about 20 bytes per non-zero plus ~0.4 GB for Python
and the libraries (an MTX folder costs more than an `.h5`, since SciPy's reader parses the text
in one go). A typical 10x sample (≤ 30 000 nuclei × ≤ 2 500 genes per nucleus) fits Lambda's
3008 MB; larger samples or a whole-atlas `.h5ad` / `.rds` need the standard runtime with more
memory, and Lambda stops a run after 15 minutes (about 150 typical samples). A Seurat file
needs roughly 80 bytes per non-zero in the reader process (1.6 GB for 20 M non-zeros), and its
temporary `.h5ad` must fit Lambda's `/tmp`. Writing the droplets is a file copy for count and
Flex outputs.

## Science and deviations from best practice

sc-prep follows the single-cell best-practices guide (Heumos et al. 2023) in keeping raw counts
untouched, reading the raw droplet matrix for ambient RNA correction, and computing the QC
covariates the guide uses (counts, genes, mitochondrial and ribosomal fractions). It deliberately
does no filtering, ambient correction or doublet detection (sc-qc's job), so that QC decisions
are made in one place with their parameters recorded. `MIN_COUNTS_PER_CELL` is only a floor
for runs that need it.

## Limitations

- It does not judge quality: empty droplets Cell Ranger called as cells, doublets and nuclei
  with ambient RNA are all kept.
- 10x Flex probe sets have no `RPS` / `RPL` probes, so `pct_counts_ribo` is 0 for Flex data
  (as for both demo datasets).
- `organism` from gene symbol case only recognises human; mouse and rat data without Ensembl
  ids, a 10x genome name or a `species` column come out `unknown`: set `ORGANISM`.
- The SPARC join needs the dataset's `sam-` / `sub-` folders or sample ids equal to the SPARC
  ids; a sample staged without them (files side by side) joins only by `sample_id`.
- Gene sets can differ between samples (e.g. the Flex probe sets of Cell Ranger 7.1 and 8.0
  give 18 082 vs 18 099 genes in Discover 475); later steps use the genes shared by all.
- An aggregated matrix split by barcode suffix gets no droplets; `.h5ad` and Seurat inputs
  never have droplets.
- Seurat: only the counts and `meta.data` are read (no reductions, graphs or images).

## Development

```sh
/opt/homebrew/bin/python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python      # pytest on the synthetic fixture, no Docker
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/
```

The fixture (`tests/make_fixture.py`) has one tiny sample per supported layout: SPARC sheets
(clashing, merged, empty and clinical columns) with Flex `sam-`/`sub-` folders and a raw-only
folder, an ARC sample with peaks, a mouse Cell Ranger count sample as MTX folders, an aggregated
MTX with two barcode suffixes, an `.h5ad` with Ensembl ids and counts in a layer, an `.h5ad`
with two samples and counts in `raw.X`, a Seurat v5 object with split layers and a bare
`dgCMatrix` (which must fail). The Seurat files in `tests/data/` were written by
`tests/data/make_rds.R` with SeuratObject. See `sample-input/README.md` for the demo data.

## Licences

Apache-2.0 for this app's code (`LICENSE`). Bundled components are BSD / MIT / Apache / MPL /
PSF licensed; nothing GPL is in the image, and R is not bundled (see `THIRD_PARTY.md`). The demo
datasets are CC BY 4.0.

## Credits

- **Scanpy** (Wolf, Angerer, Theis; scverse) computes the per-cell QC metrics.
- **anndata** (Virshup, Rybakov, Theis, Angerer, Wolf; scverse) is the file format every step
  reads and writes; sc-prep reads `.h5ad` inputs element by element with it.
- **The single-cell best-practices guide** (Heumos, Schaar, Lance et al. and the Single-cell
  Best Practices Consortium) sets the QC covariates and the division of work between prep and QC.
- **NumPy** (Harris, Millman, van der Walt, Gommers et al.), **SciPy** (Virtanen, Gommers,
  Oliphant et al. and the SciPy 1.0 contributors) and **pandas** (McKinney and the pandas
  developers) hold the arrays, sparse matrices (and read MTX files) and tables.
- **h5py** (Andrew Collette and contributors) reads 10x HDF5 matrices and writes the droplets.
- **rdata** (Ramos-Carreño, Rossi) parses Seurat `.rds` files without R.
- **Seurat / SeuratObject** (Hao, Stuart, Kowalski, Choudhary, Hoffman, Hartman, Srivastava,
  Molla, Madad, Fernandez-Granda, Satija) defines the object layout the Seurat reader follows.
- **10x Genomics Chromium / Cell Ranger** (Zheng, Terry, Belgrader et al.) defines the
  feature-barcode matrices read here.
- **openpyxl** (Eric Gazoni, Charlie Clark and contributors), hardened by **defusedxml**
  (Christian Heimes), reads the SPARC `.xlsx` sheets.
- **SPARC Data Structure** (Bandrowski, Grethe, Pilko, Gillespie, Pine, Patel, Surles-Zeigler,
  Martone) defines the `samples` / `subjects` sheets and `sam-` / `sub-` folders the metadata
  join follows.
- **Discover 475** (Arendt-Tranholm A, Sankaranarayanan I, Price TJ; CC BY 4.0), described by
  Arendt-Tranholm et al. 2026 in *Brain*, is the demo dataset all outputs were tested on.
- **Discover 480** (Sankaranarayanan I, Price TJ; CC BY 4.0) is the second demo dataset.

## Citations

Each verified on Crossref (title, full author list, year, journal, volume, issue, pages),
PubMed (PMID, by DOI search and esummary), DataCite and the Pennsieve Discover API (datasets);
full author lists in `CITATION.cff`.

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
  *Genome Biol* 2018;19:15. doi:10.1186/s13059-017-1382-0, PMID 29409532.
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: Access and store annotated data
  matrices. *J Open Source Softw* 2024;9(101):4371. doi:10.21105/joss.04371.
- Heumos L, Schaar AC, Lance C, et al. Best practices for single-cell analysis across
  modalities. *Nat Rev Genet* 2023;24(8):550–572. doi:10.1038/s41576-023-00586-w, PMID 37002403.
- Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*
  2020;585(7825):357–362. doi:10.1038/s41586-020-2649-2, PMID 32939066.
- Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific
  computing in Python. *Nat Methods* 2020;17(3):261–272. doi:10.1038/s41592-019-0686-2,
  PMID 32015543.
- McKinney W. Data Structures for Statistical Computing in Python. *Proceedings of the Python
  in Science Conference* 2010:56–61. doi:10.25080/Majora-92bf1922-00a.
- Ramos-Carreño C, Rossi T. rdata: A Python library for R datasets. *J Open Source Softw*
  2024;9(104):7540. doi:10.21105/joss.07540.
- Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and
  scalable single-cell analysis. *Nat Biotechnol* 2024;42(2):293–304.
  doi:10.1038/s41587-023-01767-y, PMID 37231261.
- Zheng GXY, Terry JM, Belgrader P, et al. Massively parallel digital transcriptional profiling
  of single cells. *Nat Commun* 2017;8:14049. doi:10.1038/ncomms14049, PMID 28091601.
- Bandrowski A, Grethe JS, Pilko A, et al. SPARC Data Structure: Rationale and Design of a FAIR
  Standard for Biomedical Research Data. bioRxiv 2021. doi:10.1101/2021.02.10.430563.
- Arendt-Tranholm A, Sankaranarayanan I, Price TJ. Single nucleus sequencing of human C2 dorsal
  root ganglia tissue from fusion surgeries (Version 1). SPARC Portal 2025.
  doi:10.26275/cxuf-ftav (Pennsieve Discover 475; CC BY 4.0).
- Arendt-Tranholm A, Sankaranarayanan I, Payne C, et al. Transcriptome of the human C2 dorsal
  root ganglia in C1–2 arthrodesis surgery: insight for neck pain. *Brain* 2026;149(4):1349–1364.
  doi:10.1093/brain/awaf370, PMID 41032656.
- Sankaranarayanan I, Price TJ. Single nucleus sequencing atlas of human dorsal root ganglion
  (Version 1). SPARC Portal 2025. doi:10.26275/9qsp-i8dh (Pennsieve Discover 480; CC BY 4.0).
- h5py 3.16.0 (BSD-3-Clause), openpyxl 3.1.5 (MIT), defusedxml 0.7.1 (PSF-2.0): software
  without a paper; see `CITATION.cff`.
