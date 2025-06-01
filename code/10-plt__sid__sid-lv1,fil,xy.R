# loading
ist <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# plotting
ps <- .plt_xy(sce, ist$clust, wcs$x)

# saving
.pdf(ps, args[[3]])
