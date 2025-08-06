# args <- list(
#     list.files("outs", "^epi", full.names=TRUE), 
#     list.files("outs", "^lv1", full.names=TRUE), 
#     "qlts/epi,lv1,krt.pdf")
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

# plotting
gs <- grep("KRT", rownames(sce), value=TRUE)
mu <- rowMeans(counts(sub[gs, ]) > 0)
ps <- lapply(gs[mu > 0.2], \(.) {
    y <- .q(logcounts(sub)[., ])
    y[!grepl("epi", sub$kid)] <- NA
    sub$foo <- y
    .plt_ps(pol, sub, "foo", lw=0) +
        theme(legend.position="none") +
        ggtitle(.) + scale_fill_gradientn(NULL,
            colors=pals::jet(), na.value="lightgrey",
            limits=c(0, 1), n.breaks=2)
})

# saving
.pdf(ps, args[[3]], 1/3)


