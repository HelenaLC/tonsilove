# args <- list(
#     list.files("outs", "^epi", full.names=TRUE),
#     list.files("outs", "^lv1", full.names=TRUE),
#     "plts/epi,lv1,krt.pdf")
args[1:2] <- lapply(args[1:2], \(.) grep("C1", ., value=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(scuttle)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
sce$kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
sce$foo <- paste(sce$sid, sce$roi, sep=".")
sub <- sce[, sce$foo == "C1.EPI03"]
sub <- logNormCounts(sub)

# aesthetics
sub$PanCK <- asinh(sub$Mean.PanCK/200)
gs <- grep("KRT", rownames(sce), value=TRUE)
ns <- gsub("KRT([0-9]+).*", "\\1", gs)
gs <- c("PanCK", gs[order(as.integer(ns))])

# plotting
ps <- lapply(gs, \(.) {
    y <- if (. == "PanCK") sub[[.]] else logcounts(sub)[., ]
    n <- sum(i <- sub$kid == "epi")
    y <- .q(ifelse(i, y, NA))
    o <- order(abs(y), na.last=FALSE)
    .plt_xy(sub[, o], unname(y)[o], split=FALSE, na=TRUE) +
        theme(legend.position="none") + ggtitle(.lab(., n)) + 
        scale_color_gradientn(
            NULL, limits=c(0, 1), n.breaks=2,
            colors=rev(pals::gnuplot()), na.value="cornsilk")
})

# saving
.pdf(ps, args[[3]], 1/2)


