# args <- list(
#     list.files("outs", "gcs", full.names=TRUE),
#     list.files("outs", "lv1", full.names=TRUE),
#     "plts/gcs,lv1,xy.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(HDF5Array)
    library(concaveman)
    library(SingleCellExperiment)
})

ps <- mapply(
    x=args[[1]], y=args[[2]], 
    SIMPLIFY=FALSE, \(x, y) {
    # loading
    sce <- readRDS(x)
    ist <- readRDS(y)
    kid <- ist$clust
    sid <- paste(sce$sid[1])
    # concave hulls
    df <- data.frame(colData(sce))
    xy <- grep("global_mm", names(df), value=TRUE)
    df <- by(df, df$roi, \(fd) concaveman(as.matrix(fd[, xy])))
    id <- rep.int(names(df), sapply(df, nrow))
    df <- data.frame(id, do.call(rbind, df))
    names(df) <- c("roi", xy <- c("x", "y"))
    mu <- df |>
        group_by(roi) |>
        summarise(across(xy, median))
    foo <- setNames(!is.na(sce$roi), colnames(sce))
    p <- .plt_xy(sce, sce$roi, sid, na=TRUE, split=FALSE) +
        geom_polygon(
            aes(x, y, group=roi), df, fill=NA, col="black",
            inherit.aes=FALSE, linewidth=0.2, linetype=2) +
        scale_color_manual(NULL, 
            values=rep(unname(pals::alphabet()), 10),
            na.value="lightgrey", breaks=levels(sce$roi)) +
        ggrepel::geom_label_repel(
            box.padding=unit(0.1, "lines"),
            label.padding=unit(0.05, "lines"),
            label.r=0.1, label.size=0.1, size=1.2,
            min.segment.length=0, segment.size=0.2,
            aes(x, y, label=roi), mu, inherit.aes=FALSE)
    p$guides$guides$colour$params$override.aes$size <- 1
    # color by distancing
    foo <- unname(sce$d)
    q <- .plt_xy(sce, foo, sid, split=FALSE)
    # add Bn layer
    s <- q$layers[[1]]$aes_params$size
    df <- q$data; df <- df[df$kid == "Bn", ]
    q <- q + ggrastr::geom_point_rast(data=df, col="cornsilk", shape=16, stroke=0, size=s)
    # add reference points
    df <- metadata(sce)$ref; names(df)[c(2, 3)] <- c("x", "y")
    q <- q + geom_point(data=df, col="black", shape=4, stroke=0.2)
    list(p, q)
}) |> Reduce(f=c)

# saving
.pdf(ps, args[[3]])
