# The Seurat test files in this folder, made with R 4.5.1 + SeuratObject 5.4.0 (Matrix 1.7):
#   docker run --rm --platform linux/amd64 -v "$PWD/tests/data:/out" rocker/r-ver:4.5.1 bash -c \
#     'Rscript -e "install.packages(\"SeuratObject\", repos=\"https://p3m.dev/cran/__linux__/noble/latest\")" && Rscript /out/make_rds.R'
# Counts are Poisson(3), set.seed(1): 8 genes x 20 cells, cells 1-10 from sample s1 and 11-20 from s2.
library(Matrix); library(SeuratObject)
set.seed(1)
genes <- c("ACTB", "CD3E", "MT-CO1", "RPL3", "SNAP25", "PRPH", "NEFH", "GAPDH")
m <- Matrix(rpois(8 * 20, 3), nrow = 8, sparse = TRUE, dimnames = list(genes, paste0("AAAC", sprintf("%04d", 1:20), "-1")))
m <- as(m, "dgCMatrix")
md <- data.frame(orig.ident = rep(c("s1", "s2"), each = 10), donor = rep(c("d1", "d2"), each = 10), row.names = colnames(m))
options(Seurat.object.assay.version = "v5")
o5 <- CreateSeuratObject(counts = m, meta.data = md)
o5[["RNA"]]$data <- log1p(o5[["RNA"]]$counts)
o5[["RNA"]] <- split(o5[["RNA"]], f = o5$orig.ident)      # layers counts.s1, counts.s2, data.s1, data.s2
saveRDS(o5, "/out/seurat_v5_split.rds")
options(Seurat.object.assay.version = "v3")
saveRDS(CreateSeuratObject(counts = m, meta.data = md), "/out/seurat_v3.rds")
saveRDS(m, "/out/dgc.rds")                                   # not a Seurat object: refused
