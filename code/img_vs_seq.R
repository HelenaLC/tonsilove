# args <- list(
#     "data/ref/sce.rds",
#     list.files("outs", "lv1", full.names=TRUE))
# 
# # dependencies
# suppressPackageStartupMessages({
#     library(dplyr)
#     library(patchwork)
#     library(SingleCellExperiment)
# })
# 
# # loading
# seq <- data.frame(colData(readRDS(args[[1]])))
# sid <- gsub(".*([A-Z][0-9]).*", "\\1", args[[2]])
# img <- lapply(args[[2]], \(.)
#     data.frame(kid=readRDS(.)$clust)) |>
#     setNames(sid) |> bind_rows(.id="sid")
# 
# # harmonize labels
# lab <- "annotation_level_1"
# kid <- as.character(seq[[lab]])
# lv2 <- "annotation_20230508"
# idx <- grep("NK|ILC", seq[[lv2]])
# kid[idx] <- "ILC"
# idx <- kid == "myeloid"
# new <- seq[[lv2]][idx]
# new[grep("DC", new)] <- "DC"
# new[grep("Mono", new)] <- "mono"
# new[grep("Neutro", new)] <- "gran"
# new[grep("Macr|Slan", new)] <- "macro"
# kid[idx] <- new
# old <- c("NBC_MBC", "CD4_T", "Cytotoxic", "Mast", "epithelial")
# new <- c("NBC/MBC", "Th", "Tc", "mast", "epi")
# idx <- match(kid, old, nomatch=0)
# kid[idx != 0] <- new[idx]
# rmv <- grepl("^pre|Cycling", kid)
# kid[rmv] <- NA; table(seq$kid <- kid)

# wrangling
df <- list(seq=seq, img=img) |>
    bind_rows(.id="typ") |>
    mutate(age=case_when(
        is.na(age) ~ gsub("[0-9]$", "", sid),
        TRUE ~ ifelse(grepl("adult", age_group), "A", "C"))) |>
    mutate(pid=case_when(
        sid == "C1" ~ "BCLL-10-T",
        sid == "A1" ~ "BCLL-15-T",
        sid == "A2" ~ "BCLL-0000",
        TRUE ~ donor_id)) |>
    mutate(pid=as.integer(factor(pid))-1) |>
    mutate(sid=case_when(!is.na(pid) ~ paste0(
        sprintf("%s-%02d-%s", age, pid, typ)), TRUE ~ sid)) |>
    mutate(kid=case_when(
        kid %in% c("PDC", "DCc", "DCp") ~ "DC",
        kid %in% c("Bl", "Bd", "GCBC") ~ "Bd/l",
        kid %in% c("Bm", "Bn", "NBC/MBC") ~ "Bn/m",
        TRUE ~ kid)) |> 
    mutate_at("kid", factor, levels=c(
        "aaa", "Bd/l", "bbb", "Bn/m", "ccc", "DC", 
        "EC", "epi", "FDC", "FRC", "gran", "ILC", 
        "macro", "mast", "mono", "PC", "Tc", "Th")) |>
    mutate(sip=case_when(typ == "seq" ~ "snRNA-seq", TRUE ~ sid))

# plotting
sub <- list(
    all=unique(df$kid),
    bcs=c("Bd/l", "Bn/m", "PC"),
    tcs=c("Tc", "Th", "ILC"),
    mye=c("DC", "gran", "mast", "mono", "macro"),
    str=c("EC", "FDC", "FRC", "epi"))

cs <- setNames(.pal, levels(df$kid))
ns <- sapply(sub, \(.) table(df$sip[df$kid %in% .]))
th <- sapply(colMaxs(ns), \(.) ceiling(./1e4)*1e4)

p <- lapply(names(sub), \(.) {
    fd <- df[df$kid %in% sub[[.]], ]
    ks <- levels(droplevels(fd$kid))
    ns <- dplyr::count(group_by(fd, sip))
    ns <- mutate(ns,
        m=case_when(n>th[.]~th[.], TRUE~n),
        l=format(n, big.mark=","))
    # cell counts
    p <- ggplot(ns, aes(sip, n)) +
        geom_col(fill="grey") +
        geom_text(
            aes(y=th[.], label=l), ns, 
            angle=90, hjust=1.1, vjust=0.5, size=1) +
        scale_y_continuous(
            limits=c(0, th[.]), breaks=c(0, th[.]),
            labels=c(0, gsub("000$", "k", th[.]))) +
        .thm_fig("bw") + theme(
            aspect.ratio=1,
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            axis.text.x=element_blank()) +
        coord_cartesian(expand=FALSE) +
        ggtitle(.lab(., nrow(fd)))
    # composition
    q <- .plt_fq(fd, "sip", "kid", hc=FALSE) +
        scale_fill_manual(values=cs, breaks=ks) +
        theme(
            aspect.ratio=2,
            plot.title=element_blank(),
            axis.title=element_blank(),
            legend.title=element_blank())
    wrap_plots(p, q, ncol=1) & scale_x_discrete(limits=\(.) rev(.))
}) |> wrap_plots(nrow=1) & theme(legend.margin=margin())

# # for each seq-based sample,
# # quantify top-5 subpopulations
# # & their relative contribution
# fqs <- df |>
#     group_by(sid, kid) |>
#     dplyr::count() |>
#     group_by(sid) |>
#     mutate(p=n/sum(n))
# top <- fqs |>
#     filter(grepl("seq", sid)) |>
#     group_by(kid) |>
#     summarise_at("p", mean) |>
#     slice_max(p, n=5) |>
#     pull(kid) |> paste()
# a <- fqs |>
#     filter(kid %in% top) |>
#     group_by(sid) |>
#     summarise_at("p", sum) |>
#     mutate(typ=df$typ[match(sid, df$sid)]) |>
#     arrange(-p) |> mutate(sid=factor(sid, sid))
# b <- group_by(a, typ) |>
#     summarise_at("p", mean) |>
#     mutate(q=round(p*100, 2)) |>
#     mutate(q=sprintf("%.2f", q))
# hl <- c("grey50", "black")[1+grepl("img", sort(unique(a$sid)))]
# q <- ggplot(a, aes(p, sid)) + 
#     geom_col(aes(fill=typ), alpha=1/3) +
#     geom_vline(aes(col=typ, xintercept=p), b, lty=2, linewidth=0.2) +
#     scale_x_continuous(limits=c(0, 1), breaks=b$p, labels=paste0(b$q, "%")) + 
#     scale_color_manual(values=c("red", "blue")) +
#     scale_fill_manual(values=c("red", "blue")) +
#     .thm_fig_d("bw") + theme(
#         aspect.ratio=2/3,
#         legend.position="none",
#         panel.grid=element_blank(),
#         axis.title=element_blank(),
#         axis.ticks=element_blank(),
#         plot.title=element_text(hjust=0.5),
#         axis.text.y=element_text(color=hl),
#         axis.text.x=element_text(color=c("red", "blue"))) +
#     coord_cartesian(expand=FALSE) + ggtitle(sprintf(
#         "contribution of 5 most frequent\nsubpopulations (%s)", 
#         paste(top, collapse=", ")))

# saving
ng <- length(gg <- list(p, q))
tf <- replicate(ng, tempfile(fileext=".pdf"), FALSE)
ws <- c(10, 5); hs <- c(6, 3)
for (. in seq_along(gg)) {
    pdf(tf[[.]], 
        width=(ws[.])/2.54, 
        height=(hs[.])/2.54)
    print(gg[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output="plts/img_vs_seq.pdf")
  