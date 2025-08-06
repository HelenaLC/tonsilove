# args <- as.list(sprintf(c(
#     "outs/gcs-%s.rds", "outs/lv1-%s.rds",
#     "outs/pol-%s.parquet", "qlts/gcs,lv1,tbm,%s.pdf"), "A2"))

# dependencies
suppressPackageStartupMessages({
    library(sp)
    library(sf)
    library(concaveman)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]]); ist <- readRDS(args[[2]])
pol <- arrow::read_parquet(args[[3]], as_data_frame=FALSE)

# wrangling
kid <- (kid <- ist$clust)[match(colnames(sce), names(kid))]
sce$kid <- factor(kid, unique(kid)[order(tolower(unique(kid)))])
ks <- c("Bd", "Bl", "Bm", "Bn", "FDC", "macro", "PC", "Th")
sce$kip <- sce$kid; sce$kip[!sce$kip %in% ks] <- NA; table(sce$kip)
xy <- grep("global_mm$", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])

# plotting
pal <- setNames(rep("white", length(ks)), ks)
pal["macro"] <- "magenta"
pal["FDC"] <- "gold"
pal["Th"] <- "blue"
pal["PC"] <- "cyan"
is <- split(colnames(sce), sce$roi)
is <- is[names(which(table(sce$roi) > 1e3))]
ps <- lapply(names(is), \(gc) {
    ch <- concaveman(xy[is[[gc]], ])
    yx <- st_coordinates(st_buffer(st_polygon(list(ch)), 0.05))
    js <- point.in.polygon(xy[,1], xy[,2], yx[,1], yx[,2])
    .plt_ps(pol, sce[, rownames(xy)[js == 1]], c="kip", id=gc) +
        scale_fill_manual(NULL, limits=names(pal),
            values=pal, na.value="lightgrey") +
        geom_polygon(
            aes(V1, V2), ch, inherit.aes=FALSE,
            fill=NA, col="black", linewidth=0.4)
})

# saving
.pdf(ps, args[[4]], 5)
