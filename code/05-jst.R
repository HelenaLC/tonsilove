# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
set.seed(250623)
pbs <- readRDS(args[[1]])
sce <- lapply(args[[2]], readRDS)
ncol(sce <- do.call(cbind, sce))

if (is.matrix(pbs)) {
    # using reference profiles
    nk <- if (wcs$sub == "tcs") 0 else seq(2, 4)
    gs <- TRUE
    nrow(pbs)
} else {
    # using feature selection only
    nk <- seq(2, 8)
    gs <- pbs
    pbs <- NULL
    length(gs)
}
# clustering
ist <- .ist(sce, gs=gs, pbs=pbs, nk=nk, ns=c(2e4, 4e4, 2e5))
ids <- intersect(letters, unique(ist$clust))
ids <- c(sort(colnames(pbs)), ids)
ist$clust <- factor(ist$clust, ids)
table(ist$clust, sce$sid)
length(unique(ist$clust))

# saving
idx <- split(colnames(sce), sce$sid)
for (sid in names(idx)) {
    jst <- lapply(ist, \(.) {
        if (!is.null(dim(.))) {
            .[intersect(idx[[sid]], rownames(.)), ] 
        } else .[intersect(idx[[sid]], names(.))]
    })
    jst$profiles <- ist$profiles
    rds <- grep(sid, args[-c(1,2)], value=TRUE)
    saveRDS(jst, rds)
}