# dependencies
set.seed(250601)
library(SingleCellExperiment)

# loading
pbs <- readRDS(args[[1]])
sce <- lapply(args[[2]], readRDS)
sce <- do.call(cbind, sce)

# analysis
ist <- .ist(sce, nk=seq(10, 15), pbs=pbs, ns=c(2e4, 4e4, 2e5))
table(ist$clust, sce$sid); length(unique(ist$clust))

# saving
idx <- split(colnames(sce), sce$sid)
for (sid in names(idx)) {
    jst <- lapply(ist, \(.) {
        if (!is.null(dim(.))) {
            .[intersect(idx[[sid]], rownames(.)), ] 
        } else .[intersect(idx[[sid]], names(.))]
    })
    jst$profiles <- ist$profiles
    rds <- grep(sid, args[-c(1, 2)], value=TRUE)
    saveRDS(jst, rds)
}