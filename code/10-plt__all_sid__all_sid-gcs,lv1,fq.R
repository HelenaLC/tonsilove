# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/gcs,lv1,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(HDF5Array)
    library(SingleCellExperiment)
})

ps <- mapply(
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
xs <- seq(-(dx <- 0.025), 1+dx, 0.05)
fd <- df |>
    filter(grepl("Bd|l", kid)) |>
    filter(!is.na(roi)) |>
    mutate(
        x=cut(d, breaks=xs),
        x=xs[as.integer(x)]+dx,
        x=factor(x, sort(unique(x))))

# quantify frequencies
fq <- fd |>
    group_by(roi, x) |>
    dplyr::count(kid) |>
    mutate(p=n/sum(n)) |>
    group_by(x, kid) |>
    summarise_at("p", mean)

# plotting
mu <- 1-mean(filter(fq, kid == "Bd")$p)
ggplot(fq, aes(x, p, fill=kid)) + 
    geom_col(alpha=2/3, position="fill", key_glyph="point") + 
    scale_fill_manual(NULL, values=c("slateblue", "lavender")) +
    scale_x_discrete("relative distance to mantle", breaks=c(0, 1)) +
    scale_y_continuous("frequency", limits=c(0, 1), n.breaks=6) +
    geom_hline(yintercept=mu, col="darkslateblue", linewidth=0.2, lty=2) +
    guides(fill=guide_legend(override.aes=list(shape=21, stroke=0, size=2))) +
    ggtitle(.lab(paste(sce$sid[1]), nrow(fd))) +
    coord_cartesian(expand=FALSE) +
    .thm_fig_c("bw") + theme(
        aspect.ratio=2/3,
        legend.position="bottom",
        panel.grid=element_blank())

})

# saving
pdf(args[[3]], onefile=TRUE, width=4/2.54, height=4/2.54)
for (p in ps) print(p); dev.off()
