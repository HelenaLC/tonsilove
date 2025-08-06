# loading
ist <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# plotting
ps <- .plt_xy(sce, ist$clust, wcs$x)
ps[-1] <- lapply(ps[-1], \(p) {
    df <- p$data
    n <- sum(df$.)
    f <- 5/log10(n)
    s <- p$layers[[1]]$aes_params$size
    p$layers[[1]]$aes_params$size <- s*f
    p
})

# saving
.pdf(ps, args[[3]])
