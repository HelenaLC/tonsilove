# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(HDF5Array)
    library(SingleCellExperiment)
}); set.seed(250626)

# loading
lys <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        sce <- readRDS(x)
        ist <- readRDS(y)
        idx <- names(kid <- ist$clust)
        kid <- kid[match(colnames(sce), idx)]
        cd <- cbind(colData(sce), kid)
        `colData<-`(sce, value=cd)
    })
pat <- ".*([A-Z][0-9]).*"
names(lys) <- gsub(pat, "\\1", names(lys))

# analysis
df <- lapply(lys, \(sce) {
    # get radial neighborhood
    cd <- colData(sce)
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(cd))
    xy <- as.matrix(cd[xy])
    nn <- nn2(xy, searchtype="radius", r=0.02, k=201)
    is <- nn$nn.idx[, -1]; is[is == 0] <- NA
    print(summary(rowSums(!is.na(is))))
    # quantify composition
    kid <- sce$kid
    id <- matrix(kid[is], nrow=ncol(sce))
    id[is.na(id)] <- ""
    names(ks) <- ks <- unique(kid)
    ns <- sapply(ks, \(k) rowSums(id == k))
    fq <- prop.table(ns, 1)
    data.frame(
        kid, fq, check.names=FALSE,
        sid=factor(sce$sid), cid=colnames(sce))
}) |> bind_rows()

# clustering
df[is.na(df)] <- 0
fd <- select(df, where(is.numeric))
km <- kmeans(fd, centers=nk <- 10)$cluster
df$ctx <- factor(km, labels=paste0("N", seq_len(nk)))
(ns <- with(df, table(ctx, sid)))

# stashing
for (sid in names(lys)) {
    idx <- match(colnames(lys[[sid]]), df$cid)
    lys[[sid]]$ctx <- df$ctx[idx]
}
(ms <- sapply(lys, \(.) table(.$ctx)))
all.equal(unname(unclass(ns)), unname(unclass(ms)))

# saving
for (. in seq_along(lys)) {
    base::saveRDS(lys[[.]], args[[2+.]])
}
