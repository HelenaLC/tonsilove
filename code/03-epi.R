# args <- list(
#     list.files("outs", "fil", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE))
# args[[3]] <- sub("fil", "epi", args[[1]])
# args <- lapply(args, \(.) grep("A2", ., value=TRUE))

# dependencies
suppressPackageStartupMessages({
    library(sp)
    library(sf)
    library(sosta)
    library(ggplot2)
    library(HDF5Array)
    library(concaveman)
    library(SpatialExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- readRDS(args[[2]])
sce$kid <- (kid <- ist$clust)[
    match(colnames(sce), names(kid))]

# identify spatial features
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
spe <- toSpatialExperiment(sce, "sid", spatialCoords=xy)
spe$foo <- "foo"; spatialCoordsNames(spe) <- c("x", "y")
est <- estimateReconstructionParametersSPE(spe, 
    marks="kid", markSelect="epi",
    imageCol="foo", plotHist=FALSE)
roi <- reconstructShapeDensityImage(spe,
    marks="kid", markSelect="epi",
    imageCol="foo", imageId="foo",
    bndw=est$bndw, thres=est$thres/10)
nrow(roi)

roi$foo <- "foo"; roi$structID <- paste0("foo", seq_len(nrow(roi))) # bug?

# exclude structures with fewer than 200 epithelial cells
idx <- assingCellsToStructures(spe, roi, imageCol="foo")
idx <- factor(idx, roi$structID)
jdx <- idx[which(spe$kid == "epi")]
rmv <- which(table(jdx) < 200)
nrow(roj <- roi[-rmv, ])

# stash assignments & shape metrics
old <- roj$structID; new <- roj$structID <- 
    sprintf("EPI%02d", seq_along(old))
spe$roi <- factor(idx, old, labels=new)
sort(table(spe$roi))

# foo <- setNames(spe$roj, colnames(spe))
# .plt_xy(spe, foo, s=0.2, split=FALSE)

# saving
base::saveRDS(spe, args[[3]])