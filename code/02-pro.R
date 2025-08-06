# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(scuttle)
    library(HDF5Array)
    library(BiocParallel)
    library(BiocSingular)
}); set.seed(250805)

# loading
sce <- readRDS(args[[1]])
pbs <- readRDS(args[[2]])

# analysis
sce <- logNormCounts(sce, BPPARAM=bp)
sel <- rownames(sce) %in% rownames(pbs)
sce <- runPCA(sce, 
    ncomponents=30, subset_row=sel,
    BSPARAM=RandomParam(), BPPARAM=bp)

# saving
rowData(sce)$sel <- sel
base::saveRDS(sce, args[[3]])