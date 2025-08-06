sce <- readRDS(args[[1]])
ps <- .plt_rgb(sce, wcs$x)
.pdf(ps, args[[2]])
