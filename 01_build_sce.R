library("DropletUtils")
library("rtracklayer")
library("lobstr")
library("here")
library("sessioninfo")

## Output directory
dir_rdata <- here("processed-data", "02_build_raw_SCE")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Locate samples
sample_paths <- dir(here("processed-data", "01_cellranger"), full.names = TRUE)

## Subset to only the "Cg" samples
sample_paths <- sample_paths[grep("-Cg-", sample_paths)]

## Code adapted from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/80e285b12b54c66126927363c725f57a1591a308/code/03_build_sce/01_build_basic_sce.R#L52-L115
## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
sce <- read10xCounts(
    samples = file.path(sample_paths, "outs", "raw_feature_bc_matrix"),
    sample.names = basename(sample_paths),
    type = "sparse",
    col.names = TRUE
)
message("RDone - ", Sys.time())

# Read 10x data and create sce - 2023-11-02 09:13:40.126057
# RDone - 2023-11-02 09:21:35.670316

## Use key similar to spe objects
sce$key <- paste0(sce$Barcode, "_", sce$Sample)

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-mm10-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

## Add sample info based on the sample IDs
sce$projection <- factor(ifelse(grepl("-neg", sce$Sample), "negative", "positive"))


## Inspect object
sce
# class: SingleCellExperiment 
# dim: 32285 3036566 
# metadata(1): Samples
# assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
#   ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(6): source type ... gene_name gene_type
# colnames(3036566): 1_AAACCCAAGAAACCCA-1 1_AAACCCAAGAAACGTC-1 ...
#   2_TTTGTTGTCTTTGCTA-1 2_TTTGTTGTCTTTGGAG-1
# colData names(4): Sample Barcode key projection
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

## Note that no empty droplets have been filtered out yet!

## Save for later
save(sce, file = file.path(dir_rdata, "Cg-sce_raw.Rdata"))

## Size in Gb
lobstr::obj_size(sce)
# 1.49 GB

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 Patched (2023-11-01 r85459)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-11-02
#  pandoc   3.1.1 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.0)
#  beachmat               2.16.0    2023-04-25 [2] Bioconductor
#  Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
#  BiocIO                 1.10.0    2023-04-25 [2] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
#  Biostrings             2.68.1    2023-05-16 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.2)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.0)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.0)
#  DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
#  DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
#  digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
#  dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
#  dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
#  DropletUtils         * 1.20.0    2023-04-25 [2] Bioconductor
#  edgeR                  3.42.4    2023-05-31 [2] Bioconductor
#  fansi                  1.0.5     2023-10-08 [2] CRAN (R 4.3.1)
#  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.0)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.4    2023-10-02 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-04-11 [2] Bioconductor
#  GenomicAlignments      1.36.0    2023-04-25 [2] Bioconductor
#  GenomicRanges        * 1.52.1    2023-10-08 [2] Bioconductor
#  ggplot2                3.4.4     2023-10-12 [2] CRAN (R 4.3.1)
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.0)
#  gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
#  HDF5Array              1.28.1    2023-05-01 [2] Bioconductor
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.0)
#  htmltools              0.5.6.1   2023-10-06 [2] CRAN (R 4.3.1)
#  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.3.0)
#  httpuv                 1.6.12    2023-10-23 [2] CRAN (R 4.3.2)
#  IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
#  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
#  later                  1.3.1     2023-05-02 [2] CRAN (R 4.3.0)
#  lattice                0.22-5    2023-10-24 [3] CRAN (R 4.3.2)
#  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.0)
#  limma                  3.56.2    2023-06-04 [2] Bioconductor
#  lobstr               * 1.1.2     2022-06-22 [2] CRAN (R 4.3.0)
#  locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.0)
#  Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.2)
#  MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [2] CRAN (R 4.3.1)
#  R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.0)
#  R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.0)
#  R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.0)
#  restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.0)
#  rhdf5                  2.44.0    2023-04-25 [2] Bioconductor
#  rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
#  Rhdf5lib               1.22.1    2023-09-10 [2] Bioconductor
#  rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.0)
#  rmote                  0.3.4     2023-05-06 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.0)
#  Rsamtools              2.16.0    2023-04-25 [2] Bioconductor
#  rtracklayer          * 1.60.1    2023-08-15 [2] Bioconductor
#  S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
#  S4Vectors            * 0.38.2    2023-09-22 [2] Bioconductor
#  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.0)
#  scuttle                1.10.3    2023-10-10 [2] Bioconductor
#  servr                  0.27      2023-05-02 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.0)
#  SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.0)
#  utf8                   1.2.4     2023-10-22 [2] CRAN (R 4.3.2)
#  vctrs                  0.6.4     2023-10-12 [2] CRAN (R 4.3.1)
#  xfun                   0.41      2023-11-01 [2] CRAN (R 4.3.2)
#  XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.0)
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
#  [1] /users/lcollado/R/4.3
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
