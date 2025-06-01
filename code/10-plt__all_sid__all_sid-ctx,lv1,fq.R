# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        ctx <- readRDS(x)
        ist <- readRDS(y)
        # wrangling
        lv1 <- ist$clust[ctx$cid]
        data.frame(ctx, lv1)
    }) |> do.call(what=rbind)

# plotting
nc <- length(unique(df$ctx))
gg <- .plt_fq(df, "ctx", "lv1", hc=TRUE) +
    scale_fill_manual(NULL, values=.pal) +
    labs(x="niche", title=NULL) +
    coord_equal(2*nc, expand=FALSE) +
    theme(
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.key.size=unit(0, "lines"))
gg$guides$guides$fill$params$override.aes$size <- 1

# saving
ggsave(args[[3]], gg, units="cm", width=4, height=4)
