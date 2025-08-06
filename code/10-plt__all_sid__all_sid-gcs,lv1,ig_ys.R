# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/gcs,lv1,ig_ys.pdf")
# args[-3] <- lapply(args[-3], \(.) grep("C1", ., value=TRUE))
# 
# # dependencies
# suppressPackageStartupMessages({
#     library(sp)
#     library(sf)
#     library(dplyr)
#     library(tidyr)
#     library(ggplot2)
#     library(scuttle)
#     library(HDF5Array)
#     library(concaveman)
#     library(SingleCellExperiment)
# })
# 
# # loading
# sce <- readRDS(args[[1]])
# ist <- readRDS(args[[2]])

# filtering
sce$kid <- ist$clust[colnames(sce)]
idx <- !is.na(sce$roi); idx <- list(
    bcs=idx & grepl("^B", sce$kid),
    tcs=idx & grepl("^Th", sce$kid))
lys <- lapply(idx, \(.) sce[, .])
sub <- do.call(cbind, lys)

# shrinking
xy <- grep("global_mm", names(colData(sub)))
xy <- as.matrix(colData(sub)[xy])
is <- split(colnames(sub), sub$roi)
cs <- lapply(is, \(.) {
    ch <- concaveman(xy[., ])
    yx <- st_coordinates(st_buffer(st_polygon(list(ch)), -0.05))
    if (nrow(yx) < 10) return(NULL)
    cs <- point.in.polygon(xy[,1], xy[,2], yx[,1], yx[,2])
    colnames(sub)[cs == 1]
})

# filtering
ex <- sapply(cs, length) < 10; ct <- cs[!ex]
roj <- setNames(rep.int(names(ct), sapply(ct, length)), unlist(ct))
sort(table(sub$roj <- roj[match(colnames(sub), names(roj))]))

# wrangling
gs <- list(
    tc=c("TGFB1", "TGFB2", "TGFB3", "IL10", "IFNG"),
    ig=c("IGHM", "IGHD", "IGHA1", "IGHG1", "IGHG2"),
    pl=c("COTL1", "MKI67", "TOP2A", "TUBB", "TUBB4B", "STMN1"))
#cpa <- sweep(counts(sub), 2, sub$Area, `/`)
cpa <- normalizeCounts(sub)
es <- t(as.matrix(cpa[unlist(gs), ]))
df <- data.frame(colData(sub), es)
fd <- df |>
    filter(!is.na(roj)) |>
    pivot_longer(all_of(unlist(unname(gs)))) |>
    mutate(name=factor(name, unlist(gs)))

# # scaling
# fd <- fd |>
#     group_by(name) |>
#     mutate_at("value", .z)

# averaging
mu <- fd |>
    group_by(roj, name) |>
    filter(grepl("^B", kid)) |>
    summarise_at("value", mean) |>
    mutate(typ=names(gs)[apply(sapply(gs, \(.) name %in% .), 1, which)])

# ordering
mx <- pivot_wider(mu)
my <- as.matrix(mx[, -1])
rownames(my) <- mx[[1]]

rnk <- mu |>
    filter(typ == "ig") |>
    mutate(value=case_when(grepl("D|M", name) ~ -value, TRUE ~ value)) |>
    group_by(name) |>
    mutate(rnk=rank(value)) |>
    group_by(roj) |>
    summarise_at("rnk", mean)
mu <- left_join(mu, rnk, by="roj")
xo <- names(sort(with(rnk, setNames(rnk, roj))))

ggplot(
    filter(mu, typ != "ig"), 
    aes(rnk, value, col=name)) +
    facet_wrap(~name, scales="free_y") +
    geom_point() + geom_smooth(formula=y~x, method="lm") 

ij <- expand.grid(gs$ig, gs$tc)
apply(ij, 1, \(.) {
    . <- as.character(.)
ggplot(
    pivot_wider(mu, id_cols="roj"), 
    aes(.data[[.[1]]], .data[[.[2]]])) +
    geom_point() + geom_smooth(formula=y~x, method="lm") 
}) |> wrap_plots()

# plotting
nx <- nrow(my)
pal <- setNames(c("purple2", "magenta2", "gold2", "cyan2", "blue2"), grep("IGH", gs, value=TRUE))
pal <- c(pal, setNames(pals::brewer.greys(6), grep("IGH", gs, invert=TRUE, value=TRUE)))
gg <- ggplot(mv, aes(roj, value, group=name, fill=name)) + 
    scale_fill_manual(NULL, values=pal) +
    # annotate("rect", alpha=0.05,
    #     fill=rep(c("gold", "blue"), nx/2)[seq(nx)],
    #     xmin=seq(0.5, nx, 1), xmax=seq(1.5, nx+1, 1),
    #     ymin=rep(-Inf, nx), ymax=rep(Inf, nx)) +
    geom_col(
        position="stack", col="white", alpha=0.8,
        linewidth=0.1, width=1, key_glyph="point") +
    geom_vline(
        xintercept=seq(0.5, 100, 1),
        col="lightgrey", linewidth=0.05) +
    geom_hline(yintercept=0, linewidth=0.1) +
    scale_y_continuous(breaks=seq(0, 10, 2), limits=c(0, 10)) +
    #scale_x_discrete(limits=.xo(my)) +
    scale_x_discrete(limits=xo) +
    coord_cartesian(expand=FALSE) +
    labs(y="mean z-scaled expr.") +
    .thm_fig_d("bw", "f") + theme(
        aspect.ratio=2/3,
        panel.grid=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))

# saving
ggsave(args[[3]], gg, units="cm", width=12, height=3.5)

foo <- rnk[match(sce$roi, names(rnk))]
p <- .plt_xy(sce, unname(foo), "C1") + 
    scale_color_gradientn(
        na.value="cornsilk", 
        colors=rev(pals::gnuplot()))
.pdf(list(p), "plts/gcs,lv1,ig_age.pdf")

pc <- data.frame(gc=rownames(my), prcomp(scale(my))$x, my)
gg <- lapply(gs, \(.)
    ggplot(pc, aes(PC1, PC2, label=gc,
        col=.z(.data[[.]]))) + ggtitle(.)) |>
    wrap_plots(nrow=1, guides="collect") & 
    geom_text(size=1) & coord_equal() &
    scale_color_gradient2(
        "z-scaled\nmean expr.", 
        low="blue", mid="grey90", high="red",
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) &
    .thm_fig_c("bw")
ggsave("plts/gcs,lv1,ig_pc.pdf", gg, unit="cm", width=15, height=3)
ggplot(pc, aes(PC1, PC2, col=IGHM, label=gc)) + geom_text() + coord_equal()
ggplot(pc, aes(PC1, PC2, col=IGHA1, label=gc)) + geom_text() + coord_equal()
ggplot(pc, aes(PC1, PC2, col=IGHG1, label=gc)) + geom_text() + coord_equal()
ggplot(pc, aes(PC1, PC2, col=IGHG2, label=gc)) + geom_text() + coord_equal()
