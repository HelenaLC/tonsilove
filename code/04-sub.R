# dependencies
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
idx <- match(
    colnames(sce), 
    names(kid <- ist$clust))
table(sce$lv1 <- kid[idx])

# define subsets
sub <- list(
    tcs=c("Th", "Tc", "ILC"),
    str=c("epi", "endo", "FDC", "FRC"),
    bcs=c("Bn", "Bm", "Bl", "Bd", "PC"),
    mye=c("macro", "mono", "mast", "gran", "DC", "PDC"))
# spot check; assure that all clusters are covered
ks <- sort(unique(kid)); ks[!ks %in% unlist(sub)]

# saving
for (. in names(sub)) {
    tmp <- sce[, sce$lv1 %in% sub[[.]]]; tmp$sub <- .
    rds <- grep(., args[-c(1,2)], value=TRUE)
    base::saveRDS(tmp, rds)
}
