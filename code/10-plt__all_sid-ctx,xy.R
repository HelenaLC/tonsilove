#args <- list(list.files("outs", "ctx", full.names=TRUE), "plts/ctx,xy.pdf")

ps <- lapply(args[[1]], \(x) {
    # loading
    sce <- readRDS(x)
    sid <- paste(sce$sid[1])
    # plotting
    ctx <- setNames(sce$ctx, colnames(sce))
    .plt_xy(sce, ctx, sid, split=FALSE)
})

# aesthetics
ps <- lapply(ps, \(.) {
    .$layers[[1]]$show.legend <- TRUE
    .$guides$guides$colour$params$override.aes$size <- 1
    . + scale_color_manual("niche", drop=FALSE, values=unname(.pal_ctx))
})

# saving
.pdf(ps, args[[2]])
