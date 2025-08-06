# loading
sce <- readRDS(args[[1]])
names(roi) <- basename(roi <- args[[2]])

# wrangling
if (length(args) == 2) {
    sce$roi <- NA
    args[[3]] <- args[[2]]
} else {
    idx <- lapply(roi, \(.) colnames(.subset_shape(sce, .align_shape(sce, .))))
    lab <- rep.int(names(roi), vapply(idx, length, numeric(1)))
    sce$roi <- lab[match(colnames(sce), unlist(idx))]
}
table(sce$roi)

# saving
base::saveRDS(sce, args[[3]])