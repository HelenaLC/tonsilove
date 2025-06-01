# dependencies
suppressPackageStartupMessages({
    library(scuttle)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
sce <- loadHDF5SummarizedExperiment(args[[1]])

# exclude cells less than 'd' from any FOV border, where
# 'd' = half the radius of an average (circular) cell
ds <- .d2b(colData(sce))
mm2 <- sce$Area.um2*1e-6
th_d <- mean(sqrt(mm2/pi))/2
ex_d <- rowAnys(ds < th_d)

# exclude low-quality cells based on
# total counts & detected features
n <- ifelse(wcs$sid == "A2", 2.5, 3.5)
.f <- \(x, n) isOutlier(x, nmads=n, log=TRUE, type="lower")
ex_n <- .f(x <- sce$nCount_RNA, n)
th_n <- attr(ex_n, "threshold")[1]
hist(log(x), n=100); abline(v=log(th_n), col="red")
ex_m <- .f(x <- sce$nFeature_RNA, n)
th_m <- attr(ex_m, "threshold")[1]
hist(log(x), n=100); abline(v=log(th_m), col="red")
ex_a <- .f(x <- sce$nCount_RNA/sce$Area.um2, n+1)
th_a <- attr(ex_a, "threshold")[1]
hist(log(x), n=100); abline(v=log(th_a), col="red")

# filtering
ex <- cbind(ex_d, ex_n, ex_m, ex_a)
round(100*colMeans(ex), 2)
sub <- sce[, !rowAnys(ex)]
round(100*ncol(sub)/ncol(sce), 2)
ncol(sce); ncol(sub)
summary(sub$Area.um2)
range(sub$nCount_RNA)
range(sub$nFeature_RNA)

# stash thresholds
ths <- list(d=th_d*1e3, n=th_n, m=th_m, a=th_a)
(metadata(sub)$ths <- ths)

# saving
base::saveRDS(sub, args[[2]])