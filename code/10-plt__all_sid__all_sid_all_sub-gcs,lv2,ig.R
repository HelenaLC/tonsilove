# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/gcs,lv2,ig_ys.pdf")
args[-3] <- lapply(args[-3], \(.) grep("C1", ., value=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(sp)
    library(sf)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggpmisc)
    library(scuttle)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- lapply(args[[2]], readRDS)

# wrangling
kid <- unlist(lapply(ist, \(.) .$clust))
sce$kid <- kid[match(colnames(sce), names(kid))]

# filtering
sub <- sce[, !is.na(sce$roi)]
sub <- logNormCounts(sub)

# shrinking
xy <- grep("global_mm", names(colData(sub)))
xy <- as.matrix(colData(sub)[xy])
is <- split(colnames(sub), sub$roi)
cs <- lapply(is, \(.) {
    ch <- concaveman(xy[., ])
    yx <- st_coordinates(st_buffer(st_polygon(list(ch)), -0.02))
    if (nrow(yx) < 10) return(NULL)
    cs <- point.in.polygon(xy[,1], xy[,2], yx[,1], yx[,2])
    colnames(sub)[cs == 1]
})

# filtering
sub$roi[-match(unlist(cs), colnames(sub))] <- NA
ex <- names(which(table(sub$roi) < 200))
sub <- sub[, !sub$roi %in% ex & !is.na(sub$roi)]

# wrangling
ks <- c(ig="^B", pl="Bd|Bgc|Bld", tc="Tfh")
gs <- list(
    ig=c("IGHM", "IGHD", "IGHA1", "IGHG1", "IGHG2"),
    pl=c("MKI67", "TOP2A", "TUBB", "COTL1", "CD40"),
    tc=c("TGFB1", "IL10", "IL17A", "IFNG", "MAF"))

# average cytokines across T cells,
# Ig's & prolif. markers across B cells
mu <- lapply(names(gs), \(.) {
    es <- logcounts(sub)[unlist(gs[.]), ]
    df <- data.frame(colData(sub), t(as.matrix(es)))
    mu <- df |>
        filter(grepl(ks[.], kid)) |>
        pivot_longer(all_of(rownames(es))) |>
        group_by(roi, name) |>
        summarise_at("value", mean)
}) |> bind_rows(.id="typ")
head(mu)

# add frequency of PCs & FDCs
fq <- data.frame(colData(sub)) |>
    group_by(roi) |>
    summarise(
        TBM=mean(kid == "macro.tbm"),
        Tfh=mean(kid == "Tfh"),
        FDC=mean(kid == "FDC"),
        PC=mean(grepl("PC", kid))) |>
    pivot_longer(-1) |>
    mutate(typ="fq")
mu <- bind_rows(mu, fq)

# get GC ranks according to -IgD/M & +IgA/G
rnk <- mu |>
    filter(name %in% gs$ig) |>
    mutate(value=case_when(
        grepl("D|M", name) ~ -value,
        TRUE ~ value)) |>
    group_by(name) |>
    mutate(rnk=rank(value)) |>
    # mutate(rnk=case_when(
    #     name == "IGHM" ~ 1*rnk,
    #     name == "IGHD" ~ 2*rnk,
    #     name == "IGHG1" ~ 3*rnk,
    #     name == "IGHA1" ~ 4*rnk,
    #     name == "IGHG2" ~ 5*rnk,
    #     TRUE ~ rnk)) |>
    group_by(roi) |>
    summarise_at("rnk", mean)
mu$rnk <- rnk$rnk[match(mu$roi, rnk$roi)]
mu$roi <- factor(mu$roi, levels(sce$roi))
mu$name <- factor(mu$name, c(unlist(gs), unique(fq$name)))

# add GC areas
as <- metadata(sub)$gcs["Area", ]
as <- data.frame(typ="fq", name="area", roi=names(as), value=as)
as$rnk <- rnk$rnk[match(as$roi, rnk$roi)]
mu <- bind_rows(mu, as)

# plotting
mv <- mu |>
    filter(name %in% gs$ig) |>
    pivot_wider()
p1 <- mapply(
    i=c("IGHD", "IGHD", "IGHD", "IGHA1", "IGHG1"),
    j=c("IGHM", "IGHA1", "IGHG1", "IGHG1", "IGHG2"),
    SIMPLIFY=FALSE, \(i, j) {
        ggplot(mv, aes(.data[[i]], .data[[j]], col=roi)) + 
            scale_color_manual(NULL, values=rep(unname(pals::polychrome()[-2]), 10)) +
            geom_point(shape=16, stroke=0, size=1) + 
            stat_poly_line(col="black", linewidth=0.2, alpha=0.2) + 
            stat_poly_eq(col="black", size=1.5, 
                label.x=Inf, label.y=Inf, vjust=1.1, hjust=1.1) 
    }) |> wrap_plots(nrow=1, guides="collect") & 
    .thm_fig_d("bw") & theme(aspect.ratio=1) &
    guides(col=guide_legend(ncol=1))

p2 <- ggplot(mu, aes(rnk, value, col=roi)) + facet_wrap(~name, scales="free_y") +
    scale_color_manual(NULL, values=rep(unname(pals::polychrome()[-2]), 10)) +
    labs(x="Ig-based maturation status", y=NULL) +
    geom_point(shape=16, stroke=0, size=1) + 
    stat_poly_line(col="black", linewidth=0.2, alpha=0.2) + 
    stat_poly_eq(col="black", size=1.5, 
        label.x=Inf, label.y=Inf, vjust=1.1, hjust=1.1) +
    .thm_fig_d("bw") + theme(aspect.ratio=1) +
    guides(col=guide_legend(ncol=1))

rnk <- distinct(mu, roi, .keep_all=TRUE)
rnk <- setNames(rnk$rnk, rnk$roi)
rnk <- rnk[match(sce$roi, names(rnk))]
p3 <- .plt_xy(sce, unname(rnk), wcs$C1, s=0.1) + 
    scale_color_gradientn(NULL, na.value="cornsilk", colors=rev(pals::gnuplot()))

# saving
pdf(args[[3]], onefile=TRUE, width=12/2.54, height=9/2.54)
for (p in list(p1, p2, p3)) print(p); dev.off()
