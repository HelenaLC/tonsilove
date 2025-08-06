# args <- list(
#     "outs/pbt.rds",
#     "data/ref/sce.rds",
#     "meta/lab/sub.json",
#     "data/ref/mty,mye.rds")

# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(scuttle)
})

# loading
pbs <- readRDS(args[[1]])
ref <- readRDS(args[[2]])

# subset myeloids
idx <- ref$annotation_level_1 %in% c("myeloid", "PDC")
table((ref <- ref[, idx])$annotation_figure_1)

# simplify annotations
sid <- "donor_id"
kid <- "annotation_20230508"
ids <- unique(ref[[kid]])
lab <- list(
    DCp=grep("PDC", ids, value=TRUE),
    macro=grep("Slan|Macro", ids, value=TRUE),
    DCc=(. <- grep("DC", ids, value=TRUE))[!grepl("PDC", .)],
    mast="Mast", mono="Monocytes", gran="Neutrophils", mye.cyc="Cycling")
idx <- match(ref[[kid]], unlist(lab))
ref$lab <- rep.int(names(lab), sapply(lab, length))[idx]
table(ref[[kid]], ref$lab)

# selection & aggregation
mgs <- findMarkers(ref, 
    groups=ref$lab, block=ref[[sid]], 
    direction="up", BPPARAM=bp)
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 100])
gs <- lapply(top, intersect, rownames(pbs))
sapply(gs, length); length(gs <- unique(unlist(gs)))
pbs <- .pbs(ref[gs, ], ids=c("lab", sid), bp=bp)

# saving
saveRDS(assay(pbs), args[[4]])
