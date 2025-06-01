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

# filtering
cs <- grep("3P$", sce$type)
sce <- logNormCounts(sce[, cs])
ref <- loadHDF5SummarizedExperiment(args[[1]])
gs <- intersect(rownames(ref), rownames(sce))
reducedDims(sce) <- list()
dim(sce <- sce[gs, ])

# saving
y <- as(assay(sce), "dgCMatrix")
assay(sce, withDimnames=FALSE) <- y
base::saveRDS(sce, args[[2]])