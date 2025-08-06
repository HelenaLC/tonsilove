# # args <- list(
# #     list.files("outs", "fil", full.names=TRUE),
# #     list.files("outs", "lv2", full.names=TRUE),
# #     "plts/fil,lv2,nn.pdf")

# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(circlize)
    library(HDF5Array)
    library(ComplexHeatmap)
    library(SingleCellExperiment)
})

# loading
.f <- \(.) split(., gsub(".*([A-Z][0-9]).*", "\\1", .))
sce <- mapply(
    x=.f(args[[1]]),
    y=.f(args[[2]]),
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        ist <- lapply(y, readRDS)
        kid <- unlist(lapply(ist, \(.) .$clust))
        kid <- kid[match(colnames(sce), names(kid))]
        sce$kid <- kid; sce
    }) |> do.call(what=cbind)

# analysis
df <- data.frame(colData(sce))
df <- by(df, df$sid, \(fd) {
    if (nrow(fd) < 200) return(NULL)
    # get radial neighborhood
    xy <- as.matrix(fd[grep("global_mm", names(fd))])
    nn <- nn2(xy, searchtype="radius", r=0.02, k=101)
    is <- nn$nn.idx[, -1]; is[is == 0] <- NA
    print(summary(rowSums(!is.na(is))))
    # quantify composition
    id <- matrix(fd$kid[is], nrow=nrow(fd))
    id[is.na(id)] <- ""
    names(ks) <- ks <- levels(factor(sce$kid))
    ns <- sapply(ks, \(k) rowSums(id == k))
    fq <- prop.table(ns, 1)
    cbind(fd[c("sid", "kid")], fq)
}) |> do.call(what=rbind)

# averaging
mu <- df |>
    group_by(kid) |>
    filter(!is.na(kid)) |>
    summarise(
        .groups="drop",
        across(where(is.numeric),
        \(.) mean(., na.rm=TRUE)))

# scaling
mv <- mu |>
    pivot_longer(where(is.numeric)) |>
    group_by(name) |> mutate_at("value", .z)

# wrangling
mx <- pivot_wider(mv)
my <- as.matrix(mx[, -1])
rownames(my) <- mx[[1]]

# plotting
pdf(args[[3]], width=8/2.54, height=8.5/2.54)
lab <- "z-scaled frequency of X in radial neighborhood of Y"
pal <- colorRamp2(c(-2.5, 0, 2.5), c("navy", "ivory", "maroon"))
set.seed(3); hm <- Heatmap(my, pal, name=lab,
    row_title=NULL, 
    column_title=NULL,
    width=unit(6, "cm"), 
    height=unit(6, "cm"),
    km=10, column_km=10,
    column_title_gp=gpar(fontsize=4),
    row_dend_gp=gpar(lwd=0.2),
    column_dend_gp=gpar(lwd=0.2),
    row_names_gp=gpar(fontsize=3),
    column_names_gp=gpar(fontsize=3),
    row_dend_width=unit(0.5, "cm"),
    column_dend_height=unit(0.5, "cm"),
    heatmap_legend_param=list(
        direction="horizontal",
        grid_height=unit(1, "mm"),
        legend_width=unit(2, "cm"),
        at=. <- seq(-2.5, 2.5, 0.5),
        labels_gp=gpar(fontsize=3),
        title_gp=gpar(fontsize=4),
        title_position="topcenter",
        labels=ifelse(. %% 2 == 0, ., "")))
draw(hm, heatmap_legend_side="top")
dev.off()
  