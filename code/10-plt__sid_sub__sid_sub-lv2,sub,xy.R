# loading
ist <- readRDS(args[[1]])
sce <- readRDS(args[[2]])

# plotting
ps <- .plt_xy(sce, ist$clust, paste(wcs$x, "-", wcs$y))
.get <- \(p) p$layers[[1]]$aes_params$size
.set <- \(p, s) { p$layers[[1]]$aes_params$size <- s; p}
ps[[1]] <- .set(ps[[1]], .get(ps[[1]])*1e5/nrow(ps[[1]]$data))
ps[-1] <- lapply(ps[-1], \(p) .set(p, .get(p)*1e4/sum(p$data$.)))

# saving
.pdf(ps, args[[3]])
