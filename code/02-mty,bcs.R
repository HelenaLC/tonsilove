# args <- list(
#     "outs/pbt.rds",
#     "data/ref/sce.rds",
#     "meta/lab/sub.json",
#     "data/ref/mty,bcs.rds")

# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(scuttle)
})

# loading
pbs <- readRDS(args[[1]])
ref <- readRDS(args[[2]])

# setup
sid <- "donor_id"
lv1 <- "annotation_level_1"
lv2 <- "annotation_20230508"

# filtering
ids <- c("PC", "preBC", "GCBC", "NBC_MBC")
ref <- ref[, ref[[lv1]] %in% ids]
table(ref[[lv1]])

# simplify annotations
lys <- list(
    "Bn"=c("NBC"),
    "Bn_IFN"=c("NBC IFN-activated"),
    "Bn.act"=c("Early GC-commited NBC", "NBC early activation", "GC-commited NBC", "preGC"),
    "Bm.cs"=c("csMBC"),
    "Bm.ncs"=c("Early MBC", "ncsMBC"),
    "Bn.p"=c("Proliferative NBC"),
    "Bm.gc"=c("Precursor MBCs", "Reactivated proliferative MBCs"),
    "Bd.np"=c("GC DZ Noproli", "DZ non proliferative"),
    "Bd.cyc"=c("DZ late G2Mphase", "DZ cell cycle exit", "DZ early Sphase", "DZ late Sphase", "DZ early G2Mphase"),
    "Bdl"=c("DZ_LZ transition"),
    "Bl"=c("LZ"),
    "Bld"=c("LZ_DZ reentry commitment", "LZ proliferative", "LZ_DZ transition"))
lab <- rep.int(names(lys), sapply(lys, length))
ref$foo <- lab[match(ref[[lv2]], unlist(lys))]
ref$foo[ref[[lv1]] == "PC"] <- "PC"; table(ref$foo)

# selection & aggregation
mgs <- findMarkers(ref, 
    groups=ref$foo, block=ref[[sid]], 
    direction="up", BPPARAM=bp)
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 100])
gs <- lapply(top, intersect, rownames(pbs))
sapply(gs, length); length(gs <- unique(unlist(gs)))
pbs <- .pbs(ref[gs, ], ids=c("foo", sid), bp=bp)

# saving
saveRDS(assay(pbs), args[[4]])
