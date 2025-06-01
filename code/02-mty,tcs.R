# args <- list(
#     "data/ref/sce.rds",
#     "data/ref/mtx.rds",
#     "data/ref/mty,tcs.rds")

# dependencies
suppressPackageStartupMessages({
    library(scran)
    library(scater)
    library(scuttle)
})

# loading
ref <- readRDS(args[[1]])
mtx <- readRDS(args[[2]])

# subset T cells
kid <- "annotation_figure_1"
idx <- c("CD4_T", "Cytotoxic")
idx <- ref$annotation_level_1 %in% idx
table((ref <- ref[, idx])[[kid]])

# simplify annotations
lys <- list(
    Thn=c("Naive", "CM Pre-non-Tfh", "CM PreTfh"),
    The=c("T-Eff-Mem", "T-Trans-Mem", "T-helper") ,
    Tfh=c("Tfh-Mem", "Tfh-LZ-GC", "Tfh T:B border", "GC-Tfh-SAP", "GC-Tfh-OX40"),
    Treg=c("Eff-Tregs", "Eff-Tregs-IL32"),
    Tfr="Tfr",
    Tcyc="cycling T",
    Tc=c("CM CD8 T","SCM CD8 T", "Naive CD8 T"),
    Trm=c("RM CD8 activated T", "RM CD8 T"),
    Tfc="CD8 Tf",
    "NK/ILC1"=c("CD16-CD56+ NK", "CD16-CD56dim NK", "CD16+CD56- NK", "ILC1"),
    ILC3=c("NKp44- ILC3","NKp44+ ILC3"))
kid <- "annotation_20230508"; sid <- "donor_id"
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
