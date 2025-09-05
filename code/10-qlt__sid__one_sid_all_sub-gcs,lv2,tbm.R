# wcs <- list(sid="A2")
# args <- as.list(sprintf(c(
#     "outs/gcs-%s.rds",
#     "outs/pol-%s.parquet",
#     "qlts/gcs,lv2,tbm,%s.pdf"), wcs$sid))
# ist <- list.files("outs", sprintf("lv2-%s", wcs$sid), full.names=TRUE)
# args <- c(args[1], list(ist), args[-1])

# dependencies
suppressPackageStartupMessages({
    library(sp)
    library(sf)
    library(arrow)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

# loading
sce <- readRDS(args[[1]])
ist <- lapply(args[[2]], readRDS)
pol <- read_parquet(args[[3]], as_data_frame=FALSE)

# wrangling
kid <- unlist(lapply(ist, \(.) .$clust))
kid <- kid[match(colnames(sce), names(kid))]
(ks <- sort(unique(sce$kid <- kid)))

ks <- c(
    pc <- grep("PC", ks, value=TRUE),
    bc <- grep("Bn|Bm", ks, value=TRUE),
    gc <- grep("Bd|Bl|gc", ks, value=TRUE),
    "FDC", "Tfh", "macro.act", "macro.tbm")
sce$kip <- sce$kid
sce$kip[!sce$kip %in% ks] <- NA
sce$kip <- factor(sce$kip)
table(sce$kip)

# aesthetics
xy <- grep("global_mm$", names(colData(sce)))
xy <- as.matrix(colData(sce)[xy])
pal <- setNames(rep("white", length(ks)), ks)
pal[bc] <- "palegreen"
pal[gc] <- "lavender"
pal["macro.act"] <- "magenta"
pal["macro.tbm"] <- "red"
pal["FDC"] <- "gold"
pal["Tfh"] <- "blue"
pal[pc] <- "cyan"

# plotting
is <- split(colnames(sce), sce$roi)
ps <- lapply(head(names(is), 10), \(gc) {
    ch <- concaveman(xy[is[[gc]], ])
    yx <- st_coordinates(st_buffer(st_polygon(list(ch)), 0.05))
    js <- point.in.polygon(xy[,1], xy[,2], yx[,1], yx[,2])
    .plt_ps(pol, sce[, rownames(xy)[js == 1]], c="kip", id=gc) +
        scale_fill_manual(NULL, limits=ks, values=pal, na.value="white") +
        geom_polygon(
            aes(V1, V2), ch, inherit.aes=FALSE,
            fill=NA, col="black", linewidth=0.2)
})

# saving
.pdf(ps, args[[4]], 5)
