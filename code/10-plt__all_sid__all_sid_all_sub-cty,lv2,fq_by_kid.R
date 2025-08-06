# args <- list(
#     list.files("outs", "cty", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/cty,lv2,fq_by_kid.pdf")

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
ps <- by(df, df$sub, \(fd) {
    id <- paste(fd$sub[1])
    .plt_fq(fd, "kid", "ctx", id) 
})
gg <- wrap_plots(ps, nrow=1, guides="collect") & 
    scale_fill_manual(NULL, values=.pal_ctx) &
    coord_equal(2*nc, expand=FALSE) & theme(
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=4),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

# saving
ggsave(args[[3]], gg, units="cm", width=15, height=3)

