#args <- list(list.files("outs", "cty", full.names=TRUE), "plts/cty,xy.pdf")

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
    . + scale_color_manual("niche", drop=FALSE, values=.pal_cty)
})

# saving
.pdf(ps, args[[2]])
