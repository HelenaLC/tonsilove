args <- list(list.files("outs", "qbs", full.names=TRUE), "plts/qbs,mg.pdf")

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

# automated selection
mgs <- lapply(lys, \(sce) {
    sce <- sce[, sce$sid == "all"]
    colnames(sce) <- names(ks) <- ks <- sce$kid
    es <- logcounts(sce)
    gs <- lapply(ks, \(i) {
        j <- setdiff(ks, i)
        fc <- es[, i]/rowMeans(es[, j])
        names(tail(sort(fc), 10))
    })
    # remove duplicates
    is <- split(seq_along(unlist(gs)), rep(names(gs), sapply(gs, length)))
    is <- lapply(is, \(.) setdiff(., which(duplicated(unlist(gs)))))
    lapply(is, \(.) unname(unlist(gs)[.]))
})

# curated selection
.mgs <- list(
    bcs=list(
        PC=c("IGHD", "IGHM", "IGHA1", "IGHG1", "IGHG2", "XBP1", "MZB1", "JCHAIN"),
        #c("JUNB", "RGS1", "RGS2", "RELA", "BTG1", "CCL3", "CCL4", "IRF4", "HMGB2", "TCL1A", "EZR", "CD38"),
        Bn=c("MS4A1", "CD19"),
        Bn_IFN=c("JUN", "FOS", "IFIT1", "IFIT3", "ISG15", "IFI44L", "MX1"),
        Bn.act=c("BIRC3", "CD40", "CD69", "CD83", "MYC", "NFKB1", "NFKBIA"),
        Bd.ncs=c(),
        Bd.cs=c(),
        Bd.np=c("CXCR5", "CXCR4"),
        Bd.p=c("MKI67", "TOP2A", "PTTG1", "STMN1", "BIRC5"),
        Bld=c("ENO1", "LDHA", "CD58"),
        Bdl=c("MAF", "IL10"),
        Bl=c("BCL2L1", "SRGN", "CD27", "CIITA", "HLA-DRB", "HLA-DRA", "HLA-DQA1", "MHC I", "CD74", "PTPRC"),
        Bm.ncs=c("TNFRSF13B"),
        Bm.cs=c("CD44", "CXCR3", "IRF4")
    ),
    tcs=list(
        c=c(),
        d=c(),
        "NK/ILC1"=c("CD3D", "CD3E", "CD3G", "NCAM1", "PRF1", "GZMB"),
        ILC3=c("NCR1", "KIT", "HLA-DPA1"),
        #ILCreg=
        a=c("GZMK", "NKG7", "CCL5"),
        Tfc=c(),
        Trm=c("CD69", "ITGAE", "ITGA1", "CD8A", "CD8B"),
        Tcn=c(),
        The=c(),
        Thn=c("TCF7", "IL7R", "CCR7", "SELL", "CD4"),
        Treg=c("FOXP3", "CTLA4", "IL2RA", "MAF", "LAG3", "IL10", "IRF4"),
        Tfh=c("ICOS", "TIGIT", "PDCD1", "CXCR5", "CXCL14"),
        Tcyc=c("MKI67", "TOP2A", "ENO1", "TUBB"),
        b=c("GPR183")
    ),
    mye=list(
        a=c("CX3CR1", "MRC1", "CSF1R"),
        mast=c("KIT", "IL1RL1", "TPSAB1/B2"),
        gran=c("FCGR3A/B", "PTGS2", "IL1B", "CXCL8", "CXCR1", "CXCR2", "CCL4/L1/L2", "S100A8", "S100A9", "IL6", "TLR2"),
        "mono_non.class"=c("CD14", "VCAN"),
        macro=c("SLC40A1", "GPNMB", "CLU", "NUPR1", "LYZ"),
        macro.act=c("CD80", "CD84", "IL18", "ITGAX"),
        TRM=c("TLR4", "CD68", "APOE", "APOC1", "CTSD", "MMP9", 
            "C1QA", "C1QB", "C1QC", "CD163", "CCL18"),
        "macro_anti.infl"=c(),
        "macro_ag.pres"=c(),
        macro.cyc=c("MKI67", "TOP2A", "CDKN1A", "AHI1", "TUBB", "PTTG1"),
        DCp=c("IL3RA", "TLR7", "CCR2", "CIITA"),
        DCc=c("HLA-DPB1", "HLA-DQB1/2", "HLA-DPA1", "CLEC10A", "LAMP3")
    ),
    str=list(
        FRC_SE=c("COL6A3", "ITGB8", "ITGA8", "PDGFRA"),
        FRC_BCZ=c("DCN", "CXCL13", "CLU", "VCAM1"),
        FDC=c("SRGN", "FN1"),
        FRC_TCZ=c("CCL19", "CXCL12", "CXCL10", "CCL2", "ANGPTL1", "COL14A1"),
        FRC_PV=c("COL5A3", "COL18A1", "NOTCH3", "ACTA2", "MYL9"),
        HEV=c("COL4A1", "CAV1", "CLEC14A", "FLT1"),
        BEC=c("CD93", "ENG", "VWF", "PECAM1", "ICAM1", "IL1R1", "IL33"),
        LEC=c("LYVE1", "PROX1", "CCL21")
    )
)
mgs$str <- list(a=unique(unname(unlist(.mgs$str))))
mgs$bcs <- .mgs$bcs
mgs$mye <- .mgs$mye
mgs$tcs <- .mgs$tcs

for (i in names(mgs)) 
    for (j in names(mgs[[i]])) 
        if (j %in% names(.mgs[[i]])) 
            mgs[[i]][[j]] <- unique(c(mgs[[i]][[j]], setdiff(.mgs[[i]][[j]], unlist(mgs))))

setdiff(unlist(mgs), rownames(lys[[1]]))
#setdiff(unlist(.mgs), rownames(lys[[1]]))
lapply(mgs, \(.) unlist(.)[duplicated(unlist(.))])
#lapply(.mgs, \(.) unlist(.)[duplicated(unlist(.))])

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
        #mutate_at("k", factor, rev(names(gs))) |>
        mutate_at("g", factor, unlist(gs))
    p <- ggplot(fd, aes(g, k, fill=y)) + aes + ggtitle(.) +
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
        #if (k %in% names(gs)) hl <- c("plain", "bold")[1 + fd$g %in% gs[[k]]]
        ggplot(fd, aes(g, s, fill=y)) + aes + ggtitle(k) +
            geom_vline(xintercept=ns+0.5, linewidth=0.2) +
            scale_y_discrete(limits=\(.) rev(.)) + guides(fill="none") +
            theme(axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5, face=hl))
    })
    c(list(p), q)
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
