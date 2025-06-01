# dependencies
suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
set.seed(250601)
pbs <- readRDS(args[[1]])
sce <- lapply(args[[2]], readRDS)
sce <- do.call(cbind, sce)

# clustering
nk <- seq(2, ifelse(is.null(pbs), 8, 5))
ist <- .ist(sce, pbs=pbs, nk=nk, ns=c(2e4, 4e4, 2e5))
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