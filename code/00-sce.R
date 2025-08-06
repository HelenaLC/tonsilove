# dependencies
suppressPackageStartupMessages({
    library(scuttle)
    library(HDF5Array)
    library(HCATonsilData)
    library(SingleCellExperiment)
})

# loading
sce <- HCATonsilData(
    assayType="RNA",
    cellType="All",
    version="2.0",
    processedCounts=FALSE)

# wrangling
cs <- grep("3P$", sce$type)
sce <- logNormCounts(sce[, cs])
reducedDims(sce) <- list()
y <- as(assay(sce), "dgCMatrix")
assay(sce, withDimnames=FALSE) <- y

# saving
base::saveRDS(sce, args[[2]])