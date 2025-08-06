# wcs <- list(sid="A1"); args <- list(
#     list.files("outs", "fil", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE))
# args[[3]] <- sub("fil", "gcs", args[[1]])
# args <- lapply(args, \(.) grep(wcs$sid, ., value=TRUE))

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

# wrangling
xy <- grep("global_mm", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
spe <- toSpatialExperiment(sce, "sid", spatialCoords=xy)
spe$foo <- "foo"; spatialCoordsNames(spe) <- c("x", "y")

# merge subpopulations of interest
kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
kid[grepl("Bd|l|FDC", kid)] <- "foo"; spe$kid <- kid

if (wcs$sid == "A2") {
    roi <- readRDS("outs/roi-A2.rds")
    roj <- roi[, grepl("GC", roi$roi)]
    is <- split(colnames(roj), roj$roi)
    yx <- lapply(is, \(.) concaveman(xy[., ]))
} else {
    # identify spatial features
    est <- estimateReconstructionParametersSPE(spe,
        marks="kid", markSelect="foo",
        imageCol="foo", plotHist=FALSE)
    roi <- reconstructShapeDensityImage(spe,
        marks="kid", markSelect="foo",
        imageCol="foo", imageId="foo",
        bndw=est$bndw/2, thres=est$thres)
    yx <- lapply(roi$sostaPolygon, as.matrix)
    yx <- lapply(yx, \(.) concaveman(., 100))
    names(yx) <- sprintf("GC%02d", seq(nrow(roi)))
}
# construct spatial feature collection
ps <- lapply(yx, \(.) st_polygon(list(.)))
nrow(roi <- st_sf(sostaPolygon=ps, id=names(yx)))

# get cell-to-structure assignments
roi$foo <- "foo"; roi$structID <- paste0("foo", seq(nrow(roi)))
idx <- assingCellsToStructures(spe, roi, imageCol="foo")
idx <- factor(idx, roi$structID)

# exclude structures with fewer than 
# 100 cells from selected subpopulations
jdx <- idx[spe$kid == "foo"]
rmv <- which(table(jdx) < 100)
nrow(roi <- roi[-rmv, ])

# stash assignments & shape metrics
old <- roi$structID; new <- roi$structID <- 
    sprintf("GC%02d", seq_along(old))
spe$roi <- factor(idx, old, labels=new)
metadata(spe)$gcs <- totalShapeMetrics(roi)
sort(table(spe$roi))

# gradients
ref <- rep(NA, nrow(roi))
names(ref) <- roi$structID
metadata(spe)$ref <- as.list(ref)
spe$d <- NA; names(spe$d) <- colnames(spe)
for (. in seq_len(nrow(roi))) {
    # get concave hulls
    ch <- concaveman(as.matrix(roi$sostaPolygon[[.]]))
    # polygon expansion
    yx <- st_coordinates(st_buffer(st_polygon(list(ch)), 0.1))
    # get 'ref'erence point (mantle midpoint)
    cs <- point.in.polygon(xy[,1], xy[,2], yx[,1], yx[,2])
    cs <- cs == 1 & spe$kid == "Bn"
    if (sum(cs) < 10) next
    ms <- colMedians(xy[cs, ])
    # get GC cell-wise distance to 'ref'
    cs <- which(spe$roi == roi$structID[.])
    ds <- apply(xy[cs, ], 1, \(.) dist(rbind(., ms)))
    metadata(spe)$ref[[roi$structID[.]]] <- ms
    spe$d[names(ds)] <- .q(ds)
}
ref <- do.call(rbind, metadata(spe)$ref)
ref <- data.frame(row.names=NULL, roi=rownames(ref), ref)
head(metadata(spe)$ref <- ref)

# saving
base::saveRDS(spe, args[[3]])
