# sample-input

sc-prep reads raw single-cell outputs. The demo datasets of the single-cell workflows are
downloaded by `scripts/fetch-sample-data.sh` in `workflow-single-cell` (Pennsieve Discover,
CC BY 4.0; requester-pays S3, so any AWS credentials work):

```sh
cd ../workflow-single-cell
AWS_PROFILE=<profile> scripts/fetch-sample-data.sh snrna      # Discover 475 + 480
```

| Folder | What | Layout sc-prep sees |
|---|---|---|
| `data/drg-c2/` | Discover 475 v1, Arendt-Tranholm, Sankaranarayanan, Price: single-nucleus RNA-seq (10x Flex, Cell Ranger 7.1 / 8.0) of C2 dorsal root ganglia from 12 neck-pain patients (doi:10.26275/cxuf-ftav) | `samples.xlsx`, `subjects.xlsx` (SPARC), `primary/sub-*/sam-*/count/sample_filtered_feature_bc_matrix.h5` + `sample_raw_feature_bc_matrix.h5`. The SPARC join gives `subject_id`, laterality, age, sex and `condition` = `subject_experimental_group` (acute / chronic). |
| `data/drg-atlas/` | Discover 480 v1, Sankaranarayanan, Price: single-nucleus atlas of human DRG, 42 donors, 83 samples (doi:10.26275/9qsp-i8dh) | As above, with derived samples one folder deeper (`sam-<donor>/sam-<donor>-L5/count/`); `condition` = `sample_experimental_group` (no pain / acute pain / unknown). |

Note: for sample `sam-UTD-CC0032-DRGR` of Discover 475 the filtered matrix is published only as
the MTX folder `count/sample_filtered_feature_bc_matrix/` (no `.h5`), which the fetch script's
`.h5` filter does not download. sc-prep reads that folder as well; fetch it to have all 12
samples.

Run the image on one of them (read-only input, outputs in `out/{cells,droplets,metadata}`):

```sh
D=../workflow-single-cell/data
docker run --rm --platform linux/amd64 --memory 3g --network none --read-only --tmpfs /tmp \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output -e MAX_CELLS_PER_SAMPLE=3000 \
  -v $D/drg-c2:/data/input:ro -v $PWD/out:/data/output sc-prep:dev
```

The outputs of these runs feed `sc-qc` (cells + droplets) and every later step (metadata); the
workflow builders keep them under `data/work/{drg-c2,drg-atlas}/prep/`.

## Your own data

Put each sample's Cell Ranger output in its own folder (`outs/` of `cellranger count` or
`cellranger-arc count`, or `outs/per_sample_outs/<sample>/` of `cellranger multi`), `.h5ad` or
Seurat `.rds` files next to them, and either the SPARC `samples.xlsx` / `subjects.xlsx` with
`sam-…` / `sub-…` folders, or a `samples.csv`:

```csv
sample_id,path,condition,donor
ctrl_1,run1/outs,control,D1
pain_1,run2/outs,pain,D2
```

An aggregated matrix (`cellranger aggr`, or a dataset's combined `filtered_feature_bc_matrix/`
whose barcodes end in `-1`, `-2`, …) is split into samples by a `barcode_suffix` column:

```csv
sample_id,path,barcode_suffix,condition
donor_a,filtered_feature_bc_matrix,1,control
donor_b,filtered_feature_bc_matrix,2,pain
```

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python    # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python        # Docker over data/fixture/linked/ -> data/output/
```
