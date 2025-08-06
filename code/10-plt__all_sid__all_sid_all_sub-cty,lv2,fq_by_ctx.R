# args <- list(
#     list.files("outs", "cty", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/cty,lv2,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
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
        data.frame(kid, sub, colData(sce)[c("ctx", "sid")])
}) |> do.call(what=rbind)
nc <- length(unique(df$ctx))

fq <- prop.table(table(df$ctx, df$sub), 1)
xs <- apply(fq, 2, \(.) rownames(fq)[. > 0.05])

# plotting
p0 <- .plt_fq(df, "ctx", "sub", "all") + 
    scale_fill_manual(NULL, values=.pal_sub)
ps <- by(df, df$sub, \(fd) {
    fq <- prop.table(table(fd$ctx))
    foo <- data.frame(table(fd$ctx, fd$kid))[[1]]
    .plt_fq(fd, "ctx", "kid", fd$sub[1], hc=FALSE, h=TRUE,
        # fade out niches with less than 10% contribution
        a=c(1/3, 1)[1+foo %in% xs[[fd$sub[1]]]]) +
        scale_fill_manual(NULL, values=.pal)
})
xs <- c("EPI", "SM", "GC", "MZ", "TCZ", "BCZ", "CTS")
ps <- lapply(c(list(p0), ps), \(p) {
    p$guides$guides$fill$params$override.aes$size <- 1
    p + scale_x_discrete(limits=xs)
})
gg <- wrap_plots(ps, nrow=1) & 
    coord_equal(2*nc, expand=FALSE) &
    theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

# saving
ggsave(args[[3]], gg, units="cm", width=15, height=3)

