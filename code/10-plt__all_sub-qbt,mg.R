#args <- list(list.files("outs", "qbt", full.names=TRUE), "plts/qbt,mg.pdf")

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
sub <- gsub(".*-([a-z]{3}).*", "\\1", args[[1]])
lys <- lapply(setNames(args[[1]], sub), readRDS)

# curated selection
mgs <- list(
    bcs=list(
        Bd.p=c(),
        Bgc.p=c("MKI67", "TOP2A", "PTTG1", "STMN1", "BIRC5"),
        Bld=c("ENO1", "LDHA", "CD58"),
        Bd.np=c("CD27", "CXCR4", "CXCR5"),
        Bm=c(),
        Bm.cs=c(),
        Bm.ncs=c("CD44", "CXCR3", "TNFRSF13B"),
        Bn_IFN=c("JUN", "FOS", "IFIT1", "IFIT3", "ISG15", "IFI44L", "MX1"),
        Bdl=c(),
        Bl=c("BCL2L1", "SRGN", "CIITA", "HLA-DRB", "HLA-DRA", "HLA-DQA1", "CD74"),
        Bn.act=c("CD40", "CD83", "BIRC3", "NFKB1", "NFKBIA", "MYC", "CD69"),
        Bn=c("MS4A1", "CD19"),
        PC_IgM=c("IGHD", "IGHM", "IGHG1", "IGHG2", "IGHA1", "XBP1", "MZB1", "IRF4", "JCHAIN"),
        PC_IgG=c(),
        PC_IgA=c()
    ),
    tcs=list(
        "NK/ILC1"=c(
            "CD8A", "CD8B", "CD4", 
            "CD3D", "CD3E", "CD3G"),
        ILC3=c(
            "NCAM1", "PRF1", "GZMB",
            "NCR1", "KIT", "HLA-DPA1"),
        Trm=c(),
        Tc=c("ITGA1", "ITGAE", "CCL5", "GZMK", "EOMES", "CXCR6", "NKG7"),
        Tcn=c(),
        Treg=c("MAF", "IRF4", "FOXP3", "CTLA4", "IL2RA", "LAG3", "IL10"),
        The=c("IL17A", "RORA", "IL18R1", "KLRB1", "CD84"),
        Thn=c("CD40LG", "TCF7", "CD69", "IL7R", "CCR7", "SELL"),
        Tfh=c("ICOS", "TIGIT", "PDCD1", "CXCR5"),
        Tcyc=c("ENO1", "MKI67", "TOP2A", "TUBB", "TUBB4B")
    ),
    mye=list(
        mast=c("KIT", "IL1RL1", "TPSAB1/B2", "RGS13"),
        gran=c("CXCL8", "CXCR1", "CXCR2", "IL1B"),
        mono.nc=c("FCGR3A/B", "S100A8", "S100A9", "VCAN"),
        mono.c=c("CD14", "IL6", "TLR2", "CD1C", "ITGAM", "NOTCH2"),
        mye.cyc=c("MKI67", "PTTG1", "TUBB", "TUBB4B", "TOP2A"),
        macro.tbm=c("CLU", "NUPR1", "APOE", "APOC1", "VCAM1"),
        macro.act=c(
            "SLC40A1", "GPNMB", "SELENOP",
            "CTSD", "MMP9", "CD68", "CD84", 
            "IL18", "CD163", "CCL18", "TLR4"),
        macro.tr=c("C1QA", "C1QB", "C1QC"),
        DC.ap=c("MRC1", "CIITA", "CLEC10A", "HLA-DPB1", "HLA-DQB1/2", "HLA-DPA1"),
        DCc=c("CD80", "LAMP3", "LY75"),
        DCp=c("IL3RA", "TLR7", "CCR2")
    ),
    str=list(
        LEC=c("LYVE1", "PROX1"),
        BEC=c("COL4A1", "CAV1", "CLEC14A", "FLT1", "CDH5"),
        BEC.act=c("CD34", "CD93", "ENG", "VWF", "PECAM1", "ICAM1", "IL1R1", "IL33"),
        FRCcts=c(
            "COL6A3", "COL1A1", "COL1A2", "COL3A1", 
            "PTGDS", "DCN", "ANGPTL1", "PDGFRA", "FOS", "CFD"),
        FRCtcz=c("COL14A1", "CCL2", "CXCL12", "CXCL10", "CCL19"),
        FRCse=c("LUM", "ITGA8", "ITGB8"),
        FRCpv=c("PDGFRB", "COL5A3", "COL18A1", "NOTCH3", "ACTA2", "MYL9"),
        FDC=c("SRGN", "CLU", "VCAM1", "CXCL13", "FN1")
    ),
    epi=list(
        epi.api1=c("KRT16", "KRT17", "KRT80", "SQSTM1", "CRYAB", "IL36G"),
        epi.api2=c("IL22RA1", "KRT23", "KRT13", "KRT4", "CEACAM6", "AREG", "IL1A"),
        epi.bas=c("KRT15", "NGFR", "ITGA3", "WNT7A"),
        epi.sup=c("KRT14", "KRT5", "STMN1", "PCNA", "WIF1", "MT1X", "ROR1", "FABP5", "WNT11"),
        epi.trn1=c("ESR1", "FAS", "JAG1", "TNFRSF10B", "PGR"),
        epi.trn2=c("CDH1", "SOX2", "HSPB1", "SFN", "LMNA"),
        epi.trn3=c("KRT7", "KRT8", "KRT18", "KRT19", "CD274", "PECAM1", "EPCAM"),
        epi.trn4=c("CCL20", "CCL22", "KRT10", "CLDN4", "REG1A", "SQLE", "MMP7")
    )
)
grep("TUBB", rownames(lys[[1]]), value=TRUE)

sce <- lys$epi
sce <- sce[, sce$sid == "all"]
colnames(sce) <- sce$kid
k <- which(colnames(sce) == "epi.trn2")
tail(sort(assay(sce)[, k]/rowMeans(assay(sce)[, -k])), 20)

lapply(sub, \(.) setdiff(
    sort(names(mgs[[.]])), 
    sort(unique(lys[[.]]$kid))))
lapply(sub, \(.) setdiff(
    sort(unique(lys[[.]]$kid)),
    sort(names(mgs[[.]]))))
setdiff(unlist(mgs), rownames(lys[[1]]))
lapply(mgs, \(.) unlist(.)[duplicated(unlist(.))])

# aesthetics
aes <- list(
    geom_tile(), 
    coord_equal(4/3, expand=FALSE),
    .thm_fig_c("minimal"), theme(
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5)),
    scale_fill_gradient2(
        "z-scaled\nmean expr.", 
        low="dodgerblue2", high="tomato2",
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)))

ps <- lapply(sub, \(.) {
    se <- lys[[.]]
    gs <- mgs[[.]]
    ns <- cumsum(sapply(gs, length))
    ns <- unique(ns); ns <- ns[-length(ns)]
    # joint
    x <- se[, se$sid == "all"]
    y <- t(logcounts(x[unlist(gs), ]))
    df <- data.frame(y, k=x$kid, check.names=FALSE)
    fd <- df |>
        pivot_longer(-k, names_to="g", values_to="y") |>
        group_by(g) |> mutate_at("y", .z) |>
        mutate_at("k", factor, rev(names(gs))) |>
        mutate_at("g", factor, unlist(gs))
    p <- ggplot(fd, aes(g, k, fill=y)) + 
        aes + ggtitle(.lab(., sum(x$ncells))) +
        geom_vline(xintercept=ns+0.5, linewidth=0.2)
    # split
    x <- se[, se$sid != "all"]
    y <- t(logcounts(x[unlist(gs), ]))
    df <- data.frame(y, s=x$sid, k=x$kid, check.names=FALSE)
    df <- mutate(df, across(where(is.numeric), .z))
    q <- by(df, df$k, \(fd) {
        k <- fd$k[1]
        fd <- fd |>
            select(-k) |>
            pivot_longer(-s, names_to="g", values_to="y") |>
            mutate_at("g", factor, unlist(gs))
        # plotting
        hl <- "plain"
        if (k %in% names(gs)) hl <- c("plain", "bold")[1 + fd$g %in% gs[[k]]]
        ggplot(fd, aes(g, s, fill=y)) + aes + ggtitle(k) +
            geom_vline(xintercept=ns+0.5, linewidth=0.2) +
            scale_y_discrete(limits=\(.) rev(.)) + guides(fill="none") +
            theme(axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5, face=hl))
    })
    list(p)#c(list(p), q)
}) |> Reduce(f=c)

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- (p <- ps[[.]])$data
    y <- ifelse(is.null(df$k), "s", "k")
    w <- length(unique(df$g))/10+max(nchar(paste(df[[y]])))/5
    h <- length(unique(df[[y]]))/10+max(nchar(paste(df$g)))/5
    pdf(tf[[.]], width=(2+w)/2.54, height=(0.5+h)/2.54)
    print(p); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[2]])
