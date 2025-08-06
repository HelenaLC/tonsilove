# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/gcs,lv1,hm.pdf")
args <- lapply(args, \(.) grep("A2", ., value=TRUE, invert=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(scater)
    library(scuttle)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
    # loading
    sce <- readRDS(x)
    ist <- readRDS(y)
    # wrangling
    sce$kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
    xy <- grep("global_mm$", names(cd <- colData(sce)))
    df <- data.frame(cd); names(df)[xy] <- c("x", "y")
    # binning
    xs <- seq(-(dx <- 0.02), 1+dx, 0.04)
    fd <- df |>
        filter(grepl("Bd|l", kid)) |>
        filter(!is.na(roi)) |>
        mutate(
            x=cut(d, breaks=xs),
            x=xs[as.integer(x)]+dx,
            x=factor(x, sort(unique(x))))
    # aggregation
    sub <- sce[, rownames(fd)]
    sub <- logNormCounts(sub)
    sub$x <- fd$x
    mu <- aggregateAcrossCells(sub, ids=sub$x, 
        use.assay.type="logcounts", statistics="mean")
    # selection
    ds <- rowDiffs(rowRanges(assay(mu)))
    gs <- names(tail(sort(ds[, 1]), 60))
    # wrangling
    gg <- data.frame(
        colData(mu), 
        t(assay(mu)[gs, ]), 
        check.names=FALSE)
    # scaling
    gg <- gg |>
        pivot_longer(all_of(gs)) |>
        group_by(name) |>
        mutate_at("value", .z)
    # ordering
    mx <- pivot_wider(select(gg, x, name, value))
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    yo <- rev(names(sort(apply(my, 2, which.max))))
    gg <- mutate_at(gg, "name", factor, yo)
    mutate(gg, n=nrow(fd))
})

# plotting
ps <- lapply(df, \(fd) {
    ggplot(fd, aes(x, name, fill=value)) + geom_tile() + 
        labs(x="relative distance to mantle", y=NULL) +
        scale_x_discrete(breaks=c(0, 1)) +
        scale_fill_gradient2(
            "z-scaled\nmean expr.", 
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2),
            low="turquoise", mid="ivory", high="purple") +
        ggtitle(.lab(paste(fd$sid[1]), fd$n[1])) +
        coord_equal(expand=FALSE) +
        .thm_fig_c("bw") + theme(
            panel.grid=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_text(size=2.5))
})

# saving
pdf(args[[3]], onefile=TRUE, width=8/2.54, height=6/2.54)
for (p in ps) print(p); dev.off()
