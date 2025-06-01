# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(scuttle)
})

# loading
ref <- readRDS(args[[1]])

# get low-level annotations
lab <- "annotation_level_1"
kid <- as.character(ref[[lab]])

# split ICLs
lv2 <- "annotation_20230508"
idx <- grep("NK|ILC", ref[[lv2]])
kid[idx] <- "ILC"

# split myeloids
idx <- kid == "myeloid"
new <- ref[[lv2]][idx]
new[grep("DC", new)] <- "DC"
new[grep("Mono", new)] <- "mono"
new[grep("Neutro", new)] <- "gran"
new[grep("Macr|Slan", new)] <- "macro"
kid[idx] <- new

# relabeling
old <- c("NBC_MBC", "CD4_T", "Cytotoxic", "Mast", "epithelial")
new <- c("NBC/MBC", "TC_CD4", "TC_CD8", "mast", "epi")
idx <- match(kid, old, nomatch=0)
kid[idx != 0] <- new[idx]

# filtering
rmv <- grepl("^pre|Cycling", kid)
ref$kid <- kid; ref <- ref[, !rmv]
dim(ref); table(ref$kid)

# aggregate counts by sample-cluster
ids <- colData(ref)[c("donor_id", "kid")]
pbs <- aggregateAcrossCells(ref, ids, BPPARAM=bp)

# library size normalization
sizeFactors(pbs) <- NULL
pbs <- logNormCounts(pbs, log=FALSE)

# average across samples
pbs <- aggregateAcrossCells(pbs, 
    ids=pbs$kid,
    statistics="mean", 
    use.assay.type="normcounts")

# saving
saveRDS(assay(pbs), args[[2]])
