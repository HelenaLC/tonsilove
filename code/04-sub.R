# dependencies
suppressPackageStartupMessages({
    library(jsonlite)
    library(SingleCellExperiment)
})

# loading
sub <- fromJSON(args[[1]])
sce <- readRDS(args[[2]])
ist <- readRDS(args[[3]])

# wrangling
idx <- match(
    colnames(sce), 
    names(kid <- ist$clust))
table(sce$lv1 <- kid[idx])

# spot check; assure that 
# all clusters are covered
ks <- sort(unique(kid))
ks[!ks %in% unlist(sub)]

# saving
for (. in names(sub)) {
    tmp <- sce[, sce$lv1 %in% sub[[.]]]
    tmp$lv1 <- factor(tmp$lv1, sub[[.]])
    rds <- grep(., args[-seq(3)], value=TRUE)
    tmp$sub <- .; base::saveRDS(tmp, rds)
}
