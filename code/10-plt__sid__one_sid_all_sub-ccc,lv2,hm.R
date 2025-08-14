# wcs <- list(x="C1"); args <- list(
#     list.files("outs", "ccc", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/ccc,lv2,hm_.pdf")
# args <- lapply(args, \(.) grep(wcs$x, ., value=TRUE, invert=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scuttle)
    library(SingleCellExperiment)
})

.sid <- \(.) gsub(".*([A-Z][0-9]).*", "\\1", .)
.sub <- \(.) gsub(".*,([a-z]{3})\\..*", "\\1", .)
.foo <- \(., by) split(., get(paste0(".", by))(.))
se <- mapply(
    x=.foo(args[[1]], "sid"),
    y=.foo(args[[2]], "sid"),
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        lys <- readRDS(x)
        ist <- lapply(.foo(y, "sub"), readRDS)
        kid <- lapply(ist, \(.) data.frame(kid=.$clust))
        kid <- bind_rows(kid, .id="sub")
        # construct 'SummarizedExperiment'
        cd <- DataFrame(kid, sid=.sid(x))
        cs <- intersect(rownames(lys[[1]]), rownames(cd))
        es <- lapply(lys, \(.)
            `rownames<-`(mx <- t(.[cs, ]),
                gsub("^(s|r)-", "", rownames(mx))))
        # characterize sender/receiver as
        # total, ligand-receptor, pathway
        rd <- data.frame(
            row.names=nm <- rownames(es[[1]]),
            typ=ifelse(grepl("total", nm), "tl",
                ifelse(grepl("-", nm), "lr", "pw")))
        se <- SummarizedExperiment(es, rowData=rd, colData=cd[cs, ])
    }) |> do.call(what=cbind)

# average by cluster
id <- colData(se)["kid"]
names(sr) <- sr <- c("s", "r")
mu <- lapply(sr, \(.) {
    #sf <- summarizeAssayByGroup(se, id, assay.type=., statistics="mean")
    ms <- summaryMarkerStats(se, groups=se$kid, assay.type=.)
    fc <- lapply(ms, \(df) {
        logFC <- log2(
            df$self.average/
                df$other.average)
        lr <- rownames(df)
        data.frame(lr, df, logFC)
    }) |> bind_rows(.id="kid")
    #fc <- lapply(ms, \(.) log2(.$self.average/.$other.average))
    #fc <- pivot_longer(bind_rows(fc, .id="kid"), -1, names_to="lr", values_to="logFC")
    # df <- data.frame(t(assay(sf)), colData(sf), check.names=FALSE)
    # fd <- pivot_longer(df, all_of(rownames(sf)), names_to="lr")
    fc |>
        #left_join(fc, by=c("kid", "lr")) |>
        mutate(typ=rowData(se)[fc$lr, "typ"]) |>
        mutate(sub=colData(se)[match(fc$kid, se$kid), "sub"])
}) |> bind_rows(.id="sr")

.hm <- \(df, ks, n1, n2, ht=FALSE, fq=0,
    val_max="logFC", val_bin="self.average") {
    df <- filter(df, kid %in% ks)
    df <- mutate(df, i=row_number())
    # exclude homotypic
    ex <- if (ht) df else 
        df |>
        group_by(lr, sr) |>
        slice_max(!!sym(val_max), n=1) |>
        group_by(lr) |>
        filter(length(unique(kid)) > 1) 
    # get top-N sender/receiver
    one <- df |>
        filter(i %in% ex$i) |>
        #group_by(kid) |>
        filter(self.detected > fq) |>
        slice_max(!!sym(val_max), n=n1, with_ties=FALSE)
    # get top-M receiver/sender
    two <- mapply(
        .sr=one$sr, .lr=one$lr, 
        SIMPLIFY=FALSE, \(.sr, .lr) {
            df |>
                filter(lr == .lr) |>
                filter(sr == setdiff(c("s", "r"), .sr)) |>
                slice_max(!!sym(val_max), n=n2, with_ties=FALSE)
        }) |> do.call(what=rbind)
    gg <- df |>
        filter(lr %in% one$lr) |>
        group_by(sr, lr) |> mutate_at(val_bin, .z) |> ungroup() |>
        mutate_at("sr", factor, c("s", "r"), c("sender", "receiver"))
    mx <- gg |> 
        filter(sr == "sender") |>
        select(all_of(c("kid", "lr", val_bin))) |>
        pivot_wider(names_from="lr", values_from=val_bin)
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    asp <- length(unique(gg$lr))/length(unique(gg$kid))
    ggplot(gg, aes(kid, lr, fill=.data[[val_bin]])) + geom_tile() +
        facet_wrap(~sr, scales="free_x") +
        scale_fill_gradient2(
            "z-scaled\nmean CCC",
            limits=c(-2.5, 2.5), n.breaks=6,
            low="cadetblue", mid="ivory", high="firebrick") +
        scale_x_discrete(limits=\(.) intersect(.xo(my), .)) +
        scale_y_discrete(limits=.yo(my)) +
        .thm_fig_c("minimal") + theme(
            aspect.ratio=asp,
            legend.position="none",
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_text(size=3),
            axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5))
}

ks <- split(mu$kid, mu$sub)
ks <- lapply(ks, unique)
ec <- grep("EC", ks$str, value=TRUE)
ilc <- grep("ILC", ks$tcs, value=TRUE)
ks <- list(
    "GC/MZ"=c(
        grep("Bn|Bd|Bl|Bgc", ks$bcs, value=TRUE), #ec,
        c("Tcyc", "Tfh", "FDC", "macro.act", "macro.tbm")),
    "CTS"=c(
        "Bn_INF", "Bm.cs", ec, ilc, "FRCpv", "FRCcts", 
        setdiff(ks$mye, c("DCc", "macro.act", "macro.tbm"))),
    "BCZ"=c(
        "FRCse", "FRCcts", 
        "Tc", "Trm", "Treg", "Tcyc", ilc,
        "PC_IgM", grep("^Bn|^Bm", ks$bcs, value=TRUE), 
        "macro.tr", "macro.act", "mye.cyc", "DC.ap", "DCp"),
    "TCZ"=c(
        "FRCtcz",
        setdiff(ks$tcs, "Tfh"),
        grep("^Bm", ks$bcs, value=TRUE),
        grep("^DC", ks$mye, value=TRUE)),
    "EPI/SM"=c(
        ks$epi, grep("PC", ks$bcs, value=TRUE),
        "FRCpv", "FRCse", ec, ilc, "Trm", "Treg",
        "gran", "mast", "macro.tr", "mye.cyc", "mono.nc")
); sapply(ks, length)

# plotting
mv <- filter(mu, typ == "lr")
ps <- lapply(names(ks), \(.) { 
    hm <- .hm(mv, ks=ks[[.]], n1=50, n2=1, ht=FALSE)
    hm + ggtitle(paste(wcs$x, ., sep="-"))
})

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    x <- ps[[.]]$data
    h <- length(unique(x$lr))/9+max(nchar(x$kid))/12
    w <- length(unique(x$kid))/5+max(nchar(x$lr))/12
    pdf(tf[[.]], 
        width=(w)/2.54, 
        height=(1+h)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
