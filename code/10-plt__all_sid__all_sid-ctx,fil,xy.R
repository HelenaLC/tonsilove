ps <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        df <- readRDS(x)
        sce <- readRDS(y)
        # plotting
        sid <- paste(sce$sid[1])
        ctx <- setNames(df$ctx, df$cid)
        .plt_xy(sce, ctx, sid, split=FALSE)
    })

# aesthetics
ps <- lapply(ps, \(.) {
    .$layers[[1]]$show.legend <- TRUE
    .$guides$guides$colour$params$override.aes$size <- 1
    . + scale_color_manual("niche", drop=FALSE, values=.pal_ctx)
})

# saving
.pdf(ps, args[[3]])
