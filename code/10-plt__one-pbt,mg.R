#args <- list("outs/pbt.rds", "plts/pbt,mg.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(scater)
    library(scuttle)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])

# restrict to selected features
gs <- list(
    epi=c("KRT5", "KRT8", "KRT14", "KRT15", "CD24"),
    EC=c("PECAM1", "CAV1", "VWF", "CD34", "CD93"),
    FRC=c("PDGFRB", "COL1A1", "COL1A2", "CCL19", "CCL21"),
    FDC=c("VCAM1", "CLU", "CXCL13", "SRGN", "FN1"),
    mast=c("KIT", "IL1RL1", "TPSAB1/B2"),
    gran=c("CXCL8", "CXCR2", "CXCR1", "PTGS2"),
    DCp=c("IL3RA", "TLR7", "CCR2"),
    DCc=c("HLA-DQB1/2", "HLA-DPA1", "CLEC10A", "ITGAX", "LAMP3"),
    mono=c("CD14", "VCAN", "FCGR3A/B", "IL1B"),
    macro=c("CD68", "TLR4", "APOE", "C1QA", "CCL18"),
    Bn=c("CD19", "MS4A1", "SELL", "CIITA", "IGHD"),
    Bm=c(),
    Bl=c("CD79A", "CD40", "CD83", "CD38"),
    Bd=c("CD27", "TOP2A", "MKI67"),
    PC=c("MZB1", "XBP1", "JCHAIN", "IGHG1", "IGHA1"),
    Th=c("CD3D", "CD3E", "CD28", "CD40LG", "KLRB1", "TCF7", "IL7R"),
    Tc=c("CD8A", "NKG7", "EOMES", "GZMK", "CXCR6", "CCL5"),
    ILC=c("KLRK1", "GNLY", "NCR1"))

unlist(gs)[duplicated(unlist(gs))]
setdiff(sce$kid, names(gs))

# wrangling
cd <- colData(sce)[c("kid", "sid")]
es <- t(logcounts(sce[unlist(gs), ]))
df <- data.frame(es, cd, check.names=FALSE)

# joint
ps <- lapply(unique(df$sid), \(.) {
    fd <- df |>
        filter(sid == .) |> 
        select(-sid) |>
        dplyr::rename(k=kid) |>
        pivot_longer(-k, names_to="g", values_to="y") |>
        group_by(g) |> mutate_at("y", .z) |>
        mutate_at("k", factor, rev(names(gs))) |>
        mutate_at("g", factor, unlist(gs))
    gg <- ggplot(fd, aes(g, k, fill=y)) + ggtitle(.) +
        geom_tile() + coord_equal(4/3, expand=FALSE) +
        scale_fill_gradient2(
            "z-scaled\nmean expr.",
            low="dodgerblue2", high="tomato2",
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
        .thm_fig_c("minimal") + theme(
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            legend.key.width=unit(0.2, "lines"),
            legend.key.height=unit(0.4, "lines"),
            axis.text.y=element_text(size=4),
            axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5))
})

# split
sub <- sce[unlist(gs), sce$sid != "all"]
es <- t(logcounts(sub)); cd <- colData(sub)
df <- data.frame(es, cd, check.names=FALSE)
df <- mutate(df, across(rownames(sub), .z))
qs <- lapply(unique(df$kid), \(.) {
    fd <- df |>
        filter(kid == .) |>
        select(-c(kid, ncells)) |>
        dplyr::rename(s=sid) |>
        pivot_longer(-s, names_to="G", values_to="y") |>
        mutate_at("G", factor, unlist(gs))
    mx <- pivot_wider(fd, names_from="s", values_from="y")
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    hl <- ifelse(rownames(sub) %in% gs[[.]], "bold", "plain")
    ggplot(fd, aes(G, s, fill=y)) + ggtitle(.) +
        geom_tile() + coord_equal(4/3, expand=FALSE) +
        scale_fill_gradient2(
            "z-scaled\nmean expr.",
            low="dodgerblue2", high="tomato2",
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) +
        scale_y_discrete(limits=.yo(my)) +
        .thm_fig_c("minimal") + theme(
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            legend.key.width=unit(0.2, "lines"),
            legend.key.height=unit(0.4, "lines"),
            axis.text.y=element_text(size=4),
            axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5, face=hl))
})

# saving
pdf(args[[2]], onefile=TRUE, width=15/2.54, height=4/2.54)
for (p in c(ps, qs)) print(p); dev.off()