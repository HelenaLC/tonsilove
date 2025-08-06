# wcs <- list(sid="C1")
# args <- as.list(sprintf(c(
#     "outs/fil-%s.rds", "outs/ccc-%s.rds", 
#     "plts/fil,ccc,xy,%s.pdf"), wcs$sid))

# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ccc <- readRDS(args[[2]])

# get nearest neighbors
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
nn <- nn2(xy, searchtype="radius", r=0.02, k=101)
summary(rowSums(nn$nn.idx[, -1] > 0))

is <- nn$nn.idx[, -1]
is[is == 0] <- NA
ys <- lapply(ccc, \(x) {
    # subset pathways
    x <- x[, -grep("-", colnames(x))]
    apply(x, 2, \(y) {
        # smooth locally
        z <- matrix(y[c(is)], nrow(is), ncol(is))
        y <- rowMeans(z, na.rm=TRUE)
        y[is.na(y)] <- 0; y
    })
    # average sender/receiver
}) |> Reduce(f=`+`)/2

# plotting
vs <- sort(colVars(ys))
yt <- tail(names(vs), 15)
ps <- lapply(yt, \(.) {
    o <- order(y <- .q(ys[, .]))
    .plt_xy(sce[, o], (y[o]), na=TRUE) + ggtitle(.) +
        scale_color_gradientn(
            "q-scaled mean\nsender/receiver",
            limits=c(0, 1), n.breaks=2,
            colors=rev(pals::gnuplot()))
})

# saving
.pdf(ps, args[[3]])
