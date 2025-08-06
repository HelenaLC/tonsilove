# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/gcs,lv1,ns.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tidytext)
    library(patchwork)
    library(SingleCellExperiment)
})

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)
        # wrangling
        sid <- paste(sce$sid[1])
        idx <- names(kid <- ist$clust)
        idx <- match(colnames(sce), idx)
        data.frame(roi=sce$roi[idx], kid=kid[idx], sid)
    }) |> 
    do.call(what=rbind) |>
    mutate(sid=factor(sid, c("C1", "A1", "A2")))

# wrangling
fd <- df |>
    filter(!is.na(roi)) |>
    filter(grepl("^B|PC|Th|FDC|macro", kid))

ns <- fd |>
    group_by(sid, roi) |>
    dplyr::count()

# plotting
.ks <- \(.) unique(.)[order(tolower(unique(.)))]
ks_sub <- .ks(fd$kid); ks_all <- .ks(df$kid)
pal <- .pal[match(ks_sub, ks_all)]
names(pal) <- ks_sub

id <- levels(df$sid)
dx <- c(C1=12e3, A1=6e3, A2=6e3)
ws <- dplyr::count(group_by(ns, sid))$n
gg <- lapply(id, \(.) {
    fd <- filter(fd, sid == .)
    fd$roi <- droplevels(fd$roi)
    p <- .plt_fq(fd, "roi", "kid", a=0.8) + ggtitle(NULL) +
        scale_fill_manual(NULL, values=pal, limits=names(pal)) +
        theme(axis.title=element_blank(), axis.text.x=element_text(size=2))
    ns <- filter(ns, sid == .)
    n <- format(nrow(fd), big.mark=",")
    m <- nlevels(fd$roi)
    q <- ggplot(ns, aes(roi, n)) + 
        geom_col(fill="grey") + 
        .thm_fig("bw") + theme(
            axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.x=element_blank()) +
        coord_cartesian(ylim=c(0, dx[.]), expand=FALSE) +
        scale_x_discrete(limits=p$scales$scales[[1]]$limits) +
        scale_y_continuous(n.breaks=4, labels=\(.) gsub("000$", "k", .)) +
        ggtitle(bquote(bold(.(.))~"(N ="~.(n)*";"~.(m)~"GCs)"))
    wrap_plots(q, p) + 
        plot_layout(ncol=1, heights=c(1, 2)) & 
        theme(plot.margin=margin(1, 3, 1, 3, unit="pt")) 
}) |> wrap_plots(nrow=1, widths=ws) + 
    plot_layout(guides="collect")

# saving
ggsave(args[[3]], gg, units="cm", width=15, height=4)
