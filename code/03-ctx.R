# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(HDF5Array)
    library(SingleCellExperiment)
}); set.seed(250530)

# neighborhoods
df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        x <- readRDS(x)
        k <- readRDS(y)$clust
        k <- k[colnames(x)]
        # get radial neighborhood
        xy <- "Center(X|Y)_global_mm"
        xy <- grep(xy, names(colData(x)))
        xy <- as.matrix(colData(x)[xy])
        nn <- nn2(xy, searchtype="radius", r=0.05, k=301)
        is <- nn$nn.idx[, -1]
        is[is == 0] <- NA
        # quantify composition
        id <- matrix(k[is], nrow=ncol(x))
        id[is.na(id)] <- ""
        names(ks) <- ks <- unique(k)
        ns <- sapply(ks, \(k) rowSums(id == k))
        fq <- prop.table(ns, 1)
        data.frame(
            sid=factor(x$sid), cid=colnames(x), 
            kid=k, fq, check.names=FALSE)
}) |> bind_rows(); df[is.na(df)] <- 0

# clustering
fd <- select(df, where(is.numeric))
km <- kmeans(fd, centers=nk <- 8)$cluster
df$ctx <- factor(km, labels=paste0("N", seq_len(nk)))
with(df, table(ctx, sid))

# saving
n <- length(s <- levels(df$sid))
for (. in seq_len(n)) {
    fd <- filter(df, sid == s[.])
    base::saveRDS(fd, args[[2+.]])
}
