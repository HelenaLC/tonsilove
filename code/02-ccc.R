# dependencies
suppressPackageStartupMessages({
    library(Matrix)
    library(scater)
    library(HDF5Array)
    library(reticulate)
    library(zellkonverter)
    library(SingleCellExperiment)
})

# setup
bin <- "~/software/mambaforge/bin/conda"
options(reticulate.conda_binary=bin)
use_condaenv("commot")
set.seed(250808)

# loading
sce <- readRDS(args[[1]])
sce <- logNormCounts(sce, BPPARAM=bp)

# wrangling
xy <- grep("global_mm$", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
reducedDim(sce, "spatial") <- xy

# run 'COMMOT'
ad <- import("anndata")
ct <- import("commot")
pd <- import("pandas")

# retrieve interactions from 'CellChatDB'
db <- ct$pp$ligand_receptor_database(
    species="human", 
    database="CellChat", 
    signaling_type=NULL)
names(db) <- c(
    "ligand", "receptor", 
    "pathway", "type")
nrow(db)

# filter for interactions
# fully represented in panel
db <- db[apply(db, 1, \(.) {
    rs <- strsplit(.["receptor"], "_")
    lr <- c(.["ligand"], unlist(rs))
    all(lr %in% rownames(sce))
}), ]
nrow(db)

# subset to features being considered
rs <- sapply(strsplit(db$receptor, "_"), .subset, 1)
nrow(sce <- sce[unique(c(db$ligand, rs)), ])

bo <- bpoptions(progressbar=TRUE)
sr <- bplapply(
    names(cs <- split(colnames(sce), sce$fov)), 
    BPPARAM=bp, BPOPTIONS=bo, \(f) {
    tryCatch(error=\(e) e, {
        # skip FOVs with very few cells
        nc <- length(cs[[f]])
        if (nc < 200) return(NULL) 
        ad <- SCE2AnnData(
            sce=sce[, cs[[f]]], 
            X_name="logcounts")
        ct$tl$spatial_communication(ad,
            dis_thr=0.02,  
            df_ligrec=db,
            heteromeric=TRUE,
            pathway_sum=TRUE,
            database_name="foo",
            heteromeric_rule="min",
            heteromeric_delimiter="_")
        list(
            ad$obsm["commot-foo-sum-sender"], 
            ad$obsm["commot-foo-sum-receiver"])
    })
})
er <- sapply(sr, inherits, "error")
ex <- sapply(sr, is.null)
table(na <- er | ex)
sr <- sr[!na]

# wrangling
i <- colnames(sce)
s <- lapply(sr, \(.) .[[1]])
r <- lapply(sr, \(.) .[[2]])
s <- do.call(rbind, unname(s))
r <- do.call(rbind, unname(r))
res <- list(s=s[i, ], r=r[i, ])
res <- lapply(res, \(.) {
    mtx <- as(as.matrix(.), "dgCMatrix")
    nms <- gsub("^(s|r)-", "", colnames(.))
    `colnames<-`(mtx, nms)
})

# saving
saveRDS(res, args[[2]])
