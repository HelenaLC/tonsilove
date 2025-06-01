args <- list("outs/pbs.rds", "plts/pbs,mg.pdf")

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
    epi_apical = c("KRT4", "KRT16", "KRT17", "KRT80", "ANXA1"),
    epi_trans1 = c("ICAM2", "CCL20", "KRT19", "ANXA2", "LYZ", "KRT8", "IDO1", "IFI6"),
    epi_trans2 = c("ICAM1", "KRT5", "PTGES", "AQP3", "WNT5A"),
    epi_basal = c("KRT14", "KRT15", "VEGFA", "EGFR", "CD24"),
    endo = c("PECAM1", "CAV1", "VWF", "CD34", "CD93"),
    FRC = c("FGFR1", "COL1A2", "PDGFRB", "ACTA2", "CCL19", "CCL21", "COL1A1", "COL3A1"),
    FDC = c("VCAM1", "CLU", "CXCL13", "SRGN", "FN1"),
    mast = c("KIT", "IL1RL1", "TPSAB1/B2"),
    gran = c("FCGR3A/B", "IL1B", "CXCL8", "CXCR2", "CXCR1", "PTGS2"),
    PDC = c("IL3RA", "TLR7", "CCR2"),
    DC = c("HLA-DPB1", "HLA-DQB1/2", "HLA-DPA1", "CLEC10A", "ITGAX", "LAMP3"),
    mono=c("CD14", "VCAN"),
    macro = c("CD14", "TLR4", "CD68", "APOE", "APOC1", "CTSD", "MMP9") ,
    TBM=c("SLC40A1", "GPNMB", "CLU", "NUPR1", "LYZ"),
    GCBC.TBM = c("SELENOP", "CCL18"),
    GCBC.LZ = c("CD79A", "CD40", "CIITA", "CD83"),
    GCBC.DZ = c("CD27", "CD38", "PTPRC", "PCNA", "TOP2A", "MKI67", "CXCR4"),
    PB = c("CD44", "CXCR3"),
    NBC = c("CD19", "MS4A1", "CD69", "SELL", "IGHD"),
    MBC = c(),
    PC_IgA1 = c("IGKC", "IGHG1", "IGHG2", "IGHA1"),
    PC_IgM = c("IRF4", "IGHM", "MZB1", "JCHAIN", "XBP1"),
    ILC = c("IL7R", "TCF7", "KLRB1", "KLRK1", "NCR1", "GNLY"),
    TC_CD4 = c("CD3D", "CD2", "CD28", "CD4", "ITGAL", "CD40LG", "TIGIT", "GATA3"),
    TC_CD8 = c("CD8A", "CXCR3", "NKG7", "EOMES", "GZMK", "CXCR6", "CCL5"))

setdiff(unlist(gs), rownames(sce))
unlist(gs)[duplicated(unlist(gs))]

# wrangling
cd <- colData(sce)[c("kid", "sid")]
es <- t(logcounts(sce[unlist(gs), ]))
df <- data.frame(es, cd, check.names=FALSE)
fd <- df |>
    filter(sid == "all") |> 
    select(-sid) |>
    dplyr::rename(k=kid) |>
    pivot_longer(-k, names_to="g", values_to="y") |>
    group_by(g) |> mutate_at("y", .z) |>
    mutate_at("k", factor, rev(names(gs))) |>
    mutate_at("g", factor, unlist(gs))

# plotting
gg <- ggplot(fd, aes(g, k, fill=y)) +
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

# saving
ggsave(gg, args[[2]], units="cm", width=15, height=4)
