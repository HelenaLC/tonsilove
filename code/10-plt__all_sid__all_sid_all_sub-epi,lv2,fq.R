# args <- list(
#     list.files("outs", "^epi", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/epi,lv2,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(HDF5Array)
    library(SingleCellExperiment)
})

.sub <- \(.) gsub(".*([a-z]{3})\\..*", "\\1", .)
.sid <- \(.) split(., gsub(".*([A-Z][0-9]).*", "\\1", .))
df <- mapply(
    x=.sid(args[[1]]),
    y=.sid(args[[2]]),
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        sid <- paste(sce$sid[1])
        ist <- lapply(y, readRDS)
        kid <- lapply(ist, \(.) .$clust)
        ns <- sapply(kid, length)
        idx <- names(kid <- unlist(kid))
        sub <- rep.int(.sub(y), ns)
        sub <- setNames(sub, idx)
        sub <- sub[match(colnames(sce), idx)]
        kid <- kid[match(colnames(sce), idx)]
        data.frame(sub, kid, colData(sce))
    }) |> do.call(what=rbind)

# wrangling
fd <- df |>
    filter(!is.na(roi)) |>
    mutate(foo=paste(sid, roi, sep=".")) |>
    mutate_at("foo", factor)

# specify regions that correspond to surface epithelium
surf <- c(
    "A1.EPI01", "A1.EPI05", "A2.EPI19", "A2.EPI23",
    "A2.EPI03", "A2.EPI04", "A2.EPI07", "A2.EPI08", 
    "A2.EPI09", "A2.EPI12", "A2.EPI13", "A2.EPI20", 
    "A2.EPI24", "A2.EPI26", "A2.EPI27")

# order by descending immune fraction
xo <- fd |>
    mutate(sup=case_when(sub %in% c("bcs", "tcs", "mye") ~ "imm", TRUE ~ sub)) |>
    group_by(foo) |> dplyr::count(sup) |> mutate(p=n/sum(n)) |>
    filter(sup == "imm") |> arrange(desc(p)) |> pull(foo)

# plotting
p0 <- .plt_fq(fd, "foo", "sub", "all", hc=TRUE, h=TRUE) + 
    scale_fill_manual(NULL, values=.pal_sub)
fq <- prop.table(table(fd$foo, fd$sub), 1)
rs <- apply(fq, 2, \(.) rownames(fq)[. > 0.05])
ps <- by(fd, fd$sub, \(x) {
    id <- paste(x$sub[1])
    p <- .plt_fq(x, "foo", "kid")
    a <- p$data$Var1 %in% rs[[id]]
    a <- ifelse(a, 1, 1/2)
    .plt_fq(x, "foo", "kid", id=id, a=1, h=TRUE) +
        scale_fill_manual(NULL, values=.pal) +
        scale_x_discrete(limits=xo)
})
ps <- lapply(c(list(p0), ps), \(p) {
    p$guides$guides$fill$params$override.aes$size <- 1
    p + scale_x_discrete(limits=xo)
})
gg <- wrap_plots(ps, nrow=1) & 
    theme(
        aspect.ratio=2,
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(
            color=.pal_sid[gsub("\\..*", "", xo)], 
            face=ifelse(xo %in% surf, "bold", "plain"),
            size=2, hjust=1, vjust=0.5))

# saving
ggsave(args[[3]], gg, units="cm", width=18, height=5)
