# dependencies
suppressPackageStartupMessages({
    library(jsonlite)
    library(SingleCellExperiment)
})

# loading
ist <- readRDS(args[[1]])
lab <- fromJSON(args[[2]])

# labelling
jst <- .jst(ist, lab)
table(ist$clust, jst$clust)

# saving
saveRDS(jst, args[[3]])