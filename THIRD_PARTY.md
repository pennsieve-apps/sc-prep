# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following
third-party components, installed from PyPI at the versions pinned in `requirements.txt`
(licences read from each package's own licence file and metadata in the built image).
All are permissively licensed (BSD, MIT, Apache-2.0, MPL-2.0, PSF); nothing GPL is in the
image. R is not in the image: Seurat `.rds` files are parsed by rdata (pure Python).

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| scanpy | 1.12.5 | QC metrics (`pp.calculate_qc_metrics`) | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata | 0.13.4 | The `.h5ad` files read (element by element) and written | BSD-3-Clause | https://github.com/scverse/anndata |
| h5py | 3.16.0 | Reads 10x HDF5 matrices, writes the droplets files (bundles HDF5, BSD-3-Clause-style) | BSD-3-Clause | https://www.h5py.org |
| rdata | 1.1.0 | Parses Seurat `.rds` files (R serialization format) without R | MIT | https://github.com/vnmabus/rdata |
| xarray | 2026.9.0 | rdata dependency (Seurat v5 cell / feature maps) | Apache-2.0 | https://xarray.dev |
| openpyxl, et_xmlfile | 3.1.5, 2.0.0 | Reads SPARC `samples.xlsx` / `subjects.xlsx` | MIT; MIT | https://foss.heptapod.net/openpyxl/openpyxl |
| defusedxml | 0.7.1 | Used by openpyxl to parse the `.xlsx` XML safely (no entity expansion) | PSF-2.0 | https://github.com/tiran/defusedxml |
| NumPy | 2.5.4 | Arrays | BSD-3-Clause (bundled parts 0BSD, MIT, Zlib, CC0-1.0) | https://numpy.org |
| SciPy | 1.18.1 | Sparse matrices (CSR / CSC), MTX reading | BSD-3-Clause | https://scipy.org |
| pandas | 3.1.0 | Tables (sample sheet, SPARC sheets, obs / var) | BSD-3-Clause | https://pandas.pydata.org |
| awslambdaric, simplejson, snapshot-restore-py | 4.2.3, 4.2.0, 1.0.0 | Lambda runtime interface client (when run as a Lambda function) | Apache-2.0; MIT OR AFL-2.1; Apache-2.0 | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| scikit-learn, joblib, threadpoolctl | 1.9.1, 1.6.0, 3.7.0 | scanpy dependencies | BSD-3-Clause | https://scikit-learn.org |
| numba, llvmlite | 0.68.1, 0.50.0 | scanpy dependencies (JIT) | BSD-2-Clause; BSD-2-Clause AND Apache-2.0 WITH LLVM-exception | https://numba.pydata.org |
| umap-learn, pynndescent | 0.5.12, 0.7.0 | scanpy dependencies (unused here) | BSD-3-Clause; BSD-2-Clause | https://github.com/lmcinnes/umap |
| statsmodels, patsy, formulaic, interface_meta, wrapt | 0.15.0, 1.0.3, 1.2.2, 2.0.1, 2.5.2 | scanpy dependencies (unused here) | BSD-3-Clause; BSD-2-Clause; MIT; MIT; BSD-2-Clause | PyPI |
| matplotlib, seaborn, contourpy, cycler, kiwisolver, fonttools, pyparsing, pillow | 3.11.2, 0.13.2, 1.4.0, 0.12.1, 1.5.1, 4.68.0, 3.3.3, 13.0.0 | scanpy dependencies (no plots are made) | Matplotlib License (PSF-based); BSD-3-Clause; BSD-3-Clause; BSD-3-Clause; BSD-3-Clause; MIT; MIT; MIT-CMU | https://matplotlib.org |
| networkx, natsort, tqdm, session-info2, legacy-api-wrap, scverse-misc, fast-array-utils | 3.7, 8.4.0, 4.70.1, 0.4.2, 1.5, 0.1.8, 1.5.1 | scanpy / anndata dependencies | BSD-3-Clause; MIT; MPL-2.0 AND MIT; MPL-2.0; MPL-2.0; BSD-3-Clause; MPL-2.0 | PyPI |
| array-api-compat, narwhals, cloudpickle, msgspec, zarr, numcodecs, donfig, google-crc32c | 1.16.0, 2.28.0, 3.1.2, 0.22.0, 3.4.2, 0.17.1, 0.8.1.post1, 1.9.0 | anndata dependencies | MIT; MIT; BSD-3-Clause; BSD-3-Clause; MIT; MIT; MIT; Apache-2.0 | PyPI |
| pydantic, pydantic-core, pydantic-settings, annotated-types, typing-inspection, python-dotenv, PyYAML | 2.14.2, 2.50.2, 2.16.2, 0.8.0, 0.4.4, 1.2.4, 6.0.3 | scverse-misc dependencies | MIT (python-dotenv BSD-3-Clause) | PyPI |
| packaging, python-dateutil, six, tzdata, typing_extensions, certifi | 26.3, 3.0.0, 1.17.0, 2026.6, 4.16.0, 2026.10.23 | dependencies | Apache-2.0 OR BSD-2-Clause; Apache-2.0/BSD dual; MIT; Apache-2.0; PSF-2.0; MPL-2.0 | PyPI |
| python:3.12.16-slim-bookworm | — | Base image (Debian 12 + CPython 3.12, pip) | PSF-2.0 (Python); Debian packages under their own free licences; pip MIT | https://hub.docker.com/_/python |

Test data in `tests/data/` (`seurat_v5_split.rds`, `seurat_v3.rds`, `dgc.rds`): synthetic Poisson
counts saved by `tests/data/make_rds.R` with R 4.5.1 and SeuratObject 5.4.0 (MIT); R and
SeuratObject are used only to make these files, not in the image.

Demo data used for testing (not in the image or the repository), both CC BY 4.0, from
Pennsieve Discover / the SPARC Portal: dataset 475, Arendt-Tranholm A, Sankaranarayanan I,
Price TJ, *Single nucleus sequencing of human C2 dorsal root ganglia tissue from fusion
surgeries* (doi:10.26275/cxuf-ftav; described by doi:10.1093/brain/awaf370); dataset 480,
Sankaranarayanan I, Price TJ, *Single nucleus sequencing atlas of human dorsal root ganglion*
(doi:10.26275/9qsp-i8dh).
