# args <- list(
#     list.files("outs", "cty", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/cty,lv1,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(HDF5Array)
    library(patchwork)
    library(SingleCellExperiment)
})

df <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- readRDS(y)
        # wrangling
        kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
        data.frame(colData(sce)[c("sid", "ctx")], kid)
    }) |> do.call(what=rbind)
nc <- length(unique(df$ctx))

# plotting
ps <- by(df, df$sid, \(fd) {
    .plt_fq(fd, "ctx", "kid", hc=TRUE) +
        ggtitle(.lab(paste(fd$sid[1]), nrow(fd))) +
        scale_fill_manual(NULL, values=.pal) +
        theme(
            axis.title=element_blank(),
            axis.text.y=element_blank(),
            legend.key.size=unit(0, "lines"))
})
xo <- ps[[1]]$scales$scales[[1]]$limits
p1 <- wrap_plots(ps, nrow=1) + 
    plot_layout(guides="collect") &
    scale_x_discrete(limits=xo)

ps <- by(df, df$ctx, \(fd) {
    .plt_fq(fd, "sid", "kid") +
        scale_x_discrete(limits=names(.pal_sid)) +
        scale_fill_manual(NULL, values=.pal) +
        labs(title=fd$ctx[1]) + theme(
            axis.title=element_blank(),
            axis.text.y=element_blank(),
            legend.key.size=unit(0, "lines"))
})
p2 <- wrap_plots(ps, nrow=1) + 
    plot_layout(guides="collect") &
    scale_x_discrete(limits=c("C1", "A1", "A2"))

a <- .plt_fq(df, "ctx", "kid", hc=TRUE) +
    scale_fill_manual("cluster", values=.pal) 
b <- .plt_fq(df, "ctx", "sid", hc=FALSE) + 
    scale_fill_manual("section", values=.pal_sid)
c <- .plt_fq(df, "sid", "ctx", hc=FALSE) + 
    scale_fill_manual("section", values=.pal_ctx)
xs <- a$scales$scales[[1]]$limits
b <- b + scale_x_discrete(limits=xs)
ps <- lapply(list(a, b, c), \(.) {
    .$guides$guides$fill$params$override.aes$size <- 1
    . + labs(title=NULL)
})
ps[[2]] <- ps[[2]] + ggtitle(.lab("all", nrow(df)))
p3 <- wrap_plots(ps, nrow=1, guides="collect") &
    theme(
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        legend.key.size=unit(0, "lines"),
        legend.spacing=unit(0, "lines"))

# saving
pdf(args[[3]], onefile=TRUE, width=10/2.54, height=3/2.54)
for (. in list(p1, p2, p3)) print(. & coord_equal(2*nc, expand=FALSE)); dev.off()
