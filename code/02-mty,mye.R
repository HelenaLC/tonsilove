# args <- list(
#     "data/ref/sce.rds",
#     "data/ref/mtx.rds",
#     "data/ref/mty,mye.rds")

# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(scuttle)
})

# loading
ref <- readRDS(args[[1]])
mtx <- readRDS(args[[2]])

# subset myeloids
idx <- ref$annotation_level_1 %in% c("myeloid", "PDC")
table((ref <- ref[, idx])$annotation_figure_1)

# simplify annotations
sid <- "donor_id"
kid <- "annotation_20230508"
ids <- unique(ref[[kid]])
lab <- list(
    PDC=grep("PDC", ids, value=TRUE),
    macro=grep("Slan|Macro", ids, value=TRUE),
    DC=(. <- grep("DC", ids, value=TRUE))[!grepl("PDC", .)],
    mast="Mast", mono="Monocytes", gran="Neutrophils", cyc="Cycling")
idx <- match(ref[[kid]], unlist(lab))
ref$lab <- rep.int(names(lab), sapply(lab, length))[idx]
table(ref[[kid]], ref$lab)

# selection & aggregation
mgs <- findMarkers(ref, groups=ref$lab, block=ref[[sid]], direction="up")
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 50])
gs <- lapply(top, intersect, rownames(mtx))
sapply(gs, length); length(gs <- unique(unlist(gs)))
pbs <- .pbs(ref[gs, ], ids=c("lab", sid))

# saving
saveRDS(assay(pbs), args[[3]])
