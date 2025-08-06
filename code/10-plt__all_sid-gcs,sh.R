#args <- list(list.files("outs", "gcs", full.names=TRUE), "plts/gcs,sh.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
df <- lapply(args[[1]], \(rds) {
    sce <- readRDS(rds)
    gcs <- metadata(sce)$gcs
    df <- data.frame(metric=rownames(gcs), gcs)
    fd <- pivot_longer(df, -1)
    data.frame(sid=sce$sid[1], fd)
}) |> do.call(what=rbind)

# wrangling
df <- df |>
    group_by(metric) |>
    mutate(val=.z(value)) |>
    mutate(sid=factor(sid, names(.pal_sid)))

# plotting
gg <- ggplot(df, aes(sid, val, fill=sid)) + 
    scale_y_continuous(limits=c(-2.5, 2.5), n.breaks=5) +
    scale_fill_manual("section", values=.pal_sid) +
    geom_hline(yintercept=0, linewidth=0.1) +
    facet_wrap(~metric, nrow=1) +
    geom_boxplot(
        alpha=2/3, key_glyph="point", 
        outlier.shape=16, outlier.size=0.1, linewidth=0.1) + 
    .thm_fig_d("bw", "f") + theme(
        aspect.ratio=3,
        legend.position="bottom",
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank()) +
    guides(fill=guide_legend(nrow=1, override.aes=list(
        alpha=1, shape=21, stroke=NA, size=1)))

# saving
ggsave(args[[2]], gg, units="cm", width=6, height=3)
