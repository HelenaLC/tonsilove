# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(arrow)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
pos <- read.csv(args[[2]])
xy <- read.csv(gzfile(args[[3]]))
yx <- filter(xy, cell %in% sce$cell)

# compute global
fs <- match(yx$fov, pos$FOV)
xs <- .mm2px(pos$X_mm[fs])
ys <- .mm2px(pos$Y_mm[fs])
dx <- max(.mm2px(pos$X_mm))
dy <- max(.mm2px(pos$Y_mm))
yx$x_global_px <- (yx$x_local_px-xs)+dx
yx$y_global_px <- dy-(yx$y_local_px+ys)

# wrangling
cd <- colData(sce)
i <- match(yx$cell, sce$cell)
j <- setdiff(names(cd), names(yx))
at <- arrow_table(cbind(yx, cd[i, j]))

# saving
write_parquet(at, args[[4]])
