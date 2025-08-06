# args <- list(
#     list.files("outs", "fil", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/fil,lv1,qc.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tidytext)
    library(HDF5Array)
    library(patchwork)
    library(SingleCellExperiment)
})

sce <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)
        # wrangling
        kid <- ist$clust
        idx <- colnames(sce)
        idx <- match(idx, names(kid))
        sce$kid <- kid[idx]
        sce
    }) |> do.call(what=cbind)

# wrangling
sub <- list(
    mye="macro|mono|mast|gran|DCc|DCp", epi="epi",
    str="EC|FDC|FRC", bcs="^B|PC", tcs="^T|ILC")
ns <- sapply(is <- lapply(sub, grep, sce$kid), length)
sce$sub[unlist(is)] <- rep.int(names(sub), ns)
table(sce$kid, sce$sub)

xs <- c("Area.um2", "cpa")
df <- data.frame(colData(sce))
df$cpa <- with(df, nCount_RNA/Area.um2)
fd <- df |>
    mutate_at("Area.um2", log10) |>
    pivot_longer(xs) |>
    mutate(name=factor(name, xs, c("log10-area(um2)", "counts/um2")))

# plotting
gg <- ggplot(fd, 
    aes(value, reorder_within(kid, value, name, median), col=sub, fill=sub)) + 
    geom_boxplot(key_glyph="point", alpha=0.2, linewidth=0.2, outlier.shape=16, outlier.size=0) + 
    facet_wrap(~name, scales="free") +
    scale_color_manual(values=.pal_sub) +
    scale_fill_manual(values=.pal_sub) +
    scale_y_reordered() +
    .thm_fig_d("bw") + theme(
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title=element_blank(),
        legend.key.size=unit(0, "pt"))

# saving
ggsave(args[[3]], gg, units="cm", width=8, height=4)
