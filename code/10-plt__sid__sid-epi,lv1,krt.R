# args <- list(
#     list.files("outs", "^epi", full.names=TRUE),
#     list.files("outs", "^lv1", full.names=TRUE),
#     "plts/epi,lv1,krt.pdf")
# args[-3] <- lapply(args[-3], \(.) grep("A2", ., value=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(scuttle)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
sce$kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
gs <- grep("KRT", rownames(sce <- logNormCounts(sce)), value=TRUE)
gs <- gs[order(as.integer(gsub("KRT([0-9]+).*", "\\1", gs)))]

# get nearest neighbors
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
nn <- nn2(xy, searchtype="radius", r=0.02, k=101)
summary(rowSums(nn$nn.idx[, -1] > 0))

# smooth locally
is <- nn$nn.idx[, -1]; is[is == 0] <- NA
es <- as.matrix(logcounts(sce)[gs, ])
ys <- apply(es, 1, \(y) {
    z <- matrix(y[c(is)], nrow(is), ncol(is))
    rowMeans(z, na.rm=TRUE)
})

# filtering
vs <- rowVars(logcounts(sce)[gs, sce$kid == "epi"])
gs <- intersect(gs, tail(names(sort(vs)), 15))

# plotting
n <- sum(i <- sce$kid == "epi")
ps <- lapply(gs, \(.) {
    y <- .q(ifelse(i, ys[, .], NA))
    o <- order(abs(y), na.last=FALSE)
    .plt_xy(sce[, o], unname(y)[o], split=FALSE, na=TRUE) +
        theme(legend.position="none") + ggtitle(.lab(., n)) +
        scale_color_gradientn(
            NULL, limits=c(0, 1), n.breaks=2,
            colors=rev(pals::gnuplot()), na.value="cornsilk")
})

# saving
.pdf(ps, args[[3]], 1/2)
