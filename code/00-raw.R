# dependencies
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(Matrix)
    library(HDF5Array)
    library(SparseArray)
    library(SingleCellExperiment)
})

# construct SCE from transcript counts & cell metadata
cd <- read.csv(gzfile(file.path(args[[1]], "metadata.csv.gz")))
mx <- readSparseCSV(file.path(args[[1]], "counts.csv.gz"), transpose=TRUE)
my <- as(mx[-1, ], "dgCMatrix"); colnames(my) <- cd$cell
sce <- SingleCellExperiment(list(counts=my), colData=cd)

# assure cell identifiers are unique across slides
colnames(sce) <- paste(did <- basename(args[[1]]), colnames(sce), sep=".")

# add local & global coordinates in mm
px <- "Center(X|Y)_(local|global)_(px)"
px <- grep(px, names(colData(sce)))
mm <- sapply(colData(sce)[px], .px2mm)
colnames(mm) <- gsub("px$", "mm", colnames(mm))
colData(sce) <- cbind(colData(sce), mm)

# move negative probe & blank code data to 'altExps'
names(gs) <- gs <- c("Negative", "SystemControl")
gs <- lapply(gs, grep, rownames(sce))
altExps(sce) <- lapply(gs, \(.) sce[., ])
sce <- sce[-unlist(gs), ]

# save data by section
md <- read.csv(args[[2]])
md <- md[md$did == did, ]
for (. in seq(nrow(md))) {
    fs <- seq(md$min[.], md$max[.])
    cs <- colnames(sce)[sce$fov %in% fs]
    sub <- sce[, cs]; sub$did <- did; sub$sid <- md$sid[.]
    saveHDF5SummarizedExperiment(sub, args[[2+.]], replace=TRUE)
}