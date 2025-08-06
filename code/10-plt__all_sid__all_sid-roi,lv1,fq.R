# args <- list(
#     list.files("outs", "roi", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/roi,lv1,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tidytext)
    library(SingleCellExperiment)
})

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)
        # wrangling
        idx <- names(kid <- ist$clust)
        idx <- match(colnames(sce), idx)
        data.frame(roi=sce$roi, kid=kid[idx])
    }) |> do.call(what=rbind)

# wrangling
fd <- df |>
    filter(!is.na(roi)) |>
    filter(grepl("GC", roi)) |>
    filter(grepl("^B|PC|Th|FDC|macro", kid))

ns <- fd |>
    group_by(roi) |>
    dplyr::count()

# plotting
ks_all <- levels(factor(df$kid))
ks_sub <- levels(factor(fd$kid))
pal <- .pal[match(ks_sub, ks_all)]
names(pal) <- ks_sub

p <- .plt_fq(fd, "roi", "kid") + theme(
    plot.title=element_blank(),
    axis.title=element_blank()) +
    scale_fill_manual(NULL, values=pal, limits=names(pal))

q <- ggplot(ns, aes(roi, n)) + 
    geom_col(fill="grey") + 
    .thm_fig() + theme(
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    coord_cartesian(ylim=c(0, 2e4), expand=FALSE) +
    scale_x_discrete(limits=p$scales$scales[[1]]$limits)

# saving
gg <- wrap_plots(q, p, ncol=1, heights=c(2, 3))
ggsave(args[[3]], gg, units="cm", width=10, height=6)
