# wcs <- list(sid="C1"); args <- c(
#     "outs/raw-%s", "outs/fil-%s.rds",
#     "outs/epi-%s.rds", "outs/gcs-%s.rds",
#     "outs/cty-%s.rds", "outs/ccc-%s.rds",
#     "outs/ist-%s.rds", "outs/lv1-%s.rds",
#     list.files(
#         "outs", full.names=TRUE,
#         paste0("(jst|lv2)-", wcs$sid)),
#     "outs/res-%s")
# args <- as.list(sprintf(args, wcs$sid))

# dependencies
suppressPackageStartupMessages({
    library(Matrix)   
    library(HDF5Array)  
    library(slingshot)  
    library(alabaster.sce)  
    library(zellkonverter)  
    library(SpatialExperiment)  
    library(SingleCellExperiment)  
})

# initialize output directory
dir <- tail(unlist(args), 1)
dir.create(dir, showWarnings=FALSE)

# load up original data
raw <- grep("raw", args, value=TRUE)
sce <- loadHDF5SummarizedExperiment(raw)
assay(sce) <- as(assay(sce), "dgCMatrix")

# quality control
fil <- readRDS(grep("fil", args, value=TRUE))
sce$fil <- colnames(sce) %in% colnames(fil)
table(sce$fil); mean(sce$fil); ncol(sce)
(metadata(sce) <- metadata(fil))

# spatial features
for (. in c("epi", "gcs")) {
    idx <- grep(paste0("^", .), basename(unlist(args)))
    roi <- readRDS(unlist(args)[idx])
    idx <- match(colnames(sce), colnames(roi))
    print(table(sce[[.]] <- roi$roi[idx], exclude=NULL))
}

# spatial contexts
ctx <- readRDS(grep("cty", args, value=TRUE))
idx <- match(colnames(sce), colnames(ctx))
table(sce$ctx <- ctx$ctx[idx])

# cell-cell communication
ccc <- readRDS(grep("ccc", args, value=TRUE))
names(ccc) <- c("sender", "receiver")
ccc <- lapply(ccc, \(mtx) {
    mty <- t(mtx[!rowAlls(as.matrix(is.na(mtx))), ])
    idx <- match(colnames(sce), colnames(mty))
    mty <- as(as.matrix(mty)[, idx], "dgCMatrix")
    `colnames<-`(mty, colnames(sce))
})
ccc <- SingleCellExperiment(ccc)
altExp(sce, "COMMOT") <- ccc

# (sub)clustering
pat <- ".*(ist|jst|lv1|lv2).*"
qat <- ".*(epi|bcs|mye|str|tcs).*"
ids <- \(.) paste(sort(unique(.)))
for (rds in grep(pat, unlist(args), value=TRUE)) {
    . <- gsub(pat, "\\1", rds)
    ist <- readRDS(rds)
    mtx <- ist$profiles
    kid <- ist$clust
    # subset
    sub <- gsub(qat, "\\1", rds)
    if (is.null(sce$sub)) sce$sub <- NA
    if (. == "lv2") {
        idx <- match(names(kid), colnames(sce))
        sce$sub[idx] <- sub
    }
    if (. %in% c("ist", "jst")) {
        # selection
        idx <- ifelse(grepl("rds", sub), "all", sub)
        sel <- rownames(sce) %in% rownames(mtx)
        rowData(sce)[[idx]] <- sel
    }
    if (is.null(sce[[.]])) {
        # new entry
        idx <- match(colnames(sce), names(kid))
        sce[[.]] <- setNames(kid[idx], colnames(sce))
    } else {
        # old entry
        old <- ids(sce[[.]])
        sce[[.]] <- setNames(paste(sce[[.]]), colnames(sce))
        sce[[.]][names(kid)] <- paste(kid)
        sce[[.]] <- factor(sce[[.]], ids(c(old, ids(kid))))
    }
}
table(sce$ist, sce$lv1); table(sce$lv2, sce$sub)
sapply(names(rowData(sce)), \(.) table(rowData(sce)[[.]]))

# saving
h5 <- file.path(dir, paste0(wcs$sid, ".h5ad"))
saveRDS(sce, file.path(dir, paste0(wcs$sid, ".rds")))
alabaster.sce::saveObject(sce, file.path(dir, wcs$sid))
zellkonverter::writeH5AD(sce, compression="gzip", h5)