# args <- list(
#     "outs/pbt.rds",
#     "data/ref/sce.rds",
#     "meta/lab/sub.json",
#     "data/ref/mty,epi.rds")

# dependencies
suppressPackageStartupMessages({
    library(jsonlite)
    library(SingleCellExperiment)
})

# loading
pbs <- readRDS(args[[1]])
sub <- fromJSON(args[[3]])

# selection
pbs <- pbs[, pbs$sid == "all"]
colnames(pbs) <- ks <- pbs$kid
a <- sub[[wcs$sub]]
b <- setdiff(ks, a)
es <- assay(pbs)
bl <- rowMeans(es[, b, drop=FALSE])
fc <- rowMaxs(sapply(a, \(.) es[, .]/bl))
sort(gs <- names(tail(sort(fc), 400)))

# saving
saveRDS(gs, args[[4]])