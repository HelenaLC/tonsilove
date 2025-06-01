# args <- list(
#     "data/ref/sce.rds",
#     "data/ref/mtx.rds",
#     "data/ref/mty,bcs.rds")

# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(scuttle)
})

# loading
ref <- readRDS(args[[1]])
mtx <- readRDS(args[[2]])

# subset B cells
sid <- "donor_id"
kid <- "annotation_20230508"
ref <- ref[, grep("BC", ref$annotation_level_1)]

# simplify annotations
lys <- list(
    "NBC"=c("NBC"),
    "NBC_IFN"=c("NBC IFN-activated"),
    "NBC_act"=c("Early GC-commited NBC", "NBC early activation", "GC-commited NBC", "preGC"),
    "MBC_cs"=c("csMBC"),
    "MBC_ncs"=c("Early MBC", "ncsMBC"),
    "NBC_p"=c("Proliferative NBC"),
    "MBC_GC"=c("Precursor MBCs", "Reactivated proliferative MBCs"),
    "DZ_np"=c("GC DZ Noproli", "DZ non proliferative"),
    "DZ_cyc"=c("DZ late G2Mphase", "DZ cell cycle exit", "DZ early Sphase", "DZ late Sphase", "DZ early G2Mphase"),
    "DZ-LZ"=c("DZ_LZ transition"),
    "LZ"=c("LZ"),
    "LZ-DZ"=c("LZ_DZ reentry commitment", "LZ proliferative", "LZ_DZ transition"))
lab <- rep.int(names(lys), sapply(lys, length))
idx <- match(ref[[kid]], unlist(lys))
table(ref[[kid <- "lab"]] <- lab[idx])

# selection & aggregation
mgs <- findMarkers(ref, groups=ref$lab, block=ref[[sid]], direction="up")
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 50])
gs <- lapply(top, intersect, rownames(mtx))
sapply(gs, length); length(gs <- unique(unlist(gs)))
pbs <- .pbs(ref[gs, ], ids=c("lab", sid))

# saving
saveRDS(assay(pbs), args[[3]])
