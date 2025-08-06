# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/gcs,lv1,ig_xy.pdf")
args[-3] <- lapply(args[-3], \(.) grep("C1", ., value=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scuttle)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])

# wrangling
sce <- logNormCounts(sce)
sce$kid <- ist$clust[colnames(sce)]

# concave hulls
df <- filter(data.frame(colData(sce)), !is.na(roi))
xy <- grep("global_mm", names(df), value=TRUE)
df <- by(df, df$roi, \(fd) concaveman(as.matrix(fd[, xy])))
id <- rep.int(names(df), sapply(df, nrow))
df <- data.frame(id, do.call(rbind, df))
names(df) <- c("roi", "x", "y")
ch <- geom_polygon(
    aes(x, y, group=roi), df, fill=NA, col="red",
    inherit.aes=FALSE, linewidth=0.2, linetype=2)

# plotting
is <- grepl("^B|PC", sce$kid)
js <- is & !is.na(sce$roi)
gs <- c("IGHM", "IGHD", "IGHA1", "IGHG1", "IGHG2")
ps <- lapply(gs, \(x) {
    lapply(list(is, js), \(y) {
        es <- logcounts(sce)[x, ]
        es <- .q(ifelse(y, es, NA))
        o <- order(es, na.last=FALSE)
        pal <- scale_color_gradientn(x, n.breaks=2, na.value="cornsilk", colors=rev(pals::gnuplot()))
        .plt_xy(sce[, o], unname(es)[o], split=FALSE) + ch + pal + ggtitle(NULL)
    })
}) |> Reduce(f=c)
.pdf(ps, args[[3]])
