#args <- list(list.files("outs", "^epi", full.names=TRUE), "plts/epi,xy.pdf")

# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(ggrepel)
    library(HDF5Array)
    library(concaveman)
    library(SpatialExperiment)
    library(SingleCellExperiment)
})

ps <- lapply(args[[1]], \(x) {
    # loading
    sce <- readRDS(x)
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
    # plotting
    foo <- setNames(sce$roi, colnames(sce))
    p <- .plt_xy(sce, foo, sid, na=TRUE, split=FALSE) + 
        scale_color_manual(
            NULL, na.value="lightgrey", 
            breaks=\(.) setdiff(., NA), 
            values=unname(pals::polychrome()[-2])) +
        geom_polygon(
            aes(x, y, group=roi), df, fill=NA, col="black",
            inherit.aes=FALSE, linewidth=0.2, linetype=2) +
        ggrepel::geom_label_repel(
            box.padding=unit(0.2, "lines"),
            label.padding=unit(0.1, "lines"),
            label.r=0.1, label.size=0.1, size=2,
            min.segment.length=0, segment.size=0.2,
            aes(x, y, label=roi), mu, inherit.aes=FALSE)
    p$guides$guides$colour$params$override.aes$size <- 1
    p
})

# saving
.pdf(ps, args[[2]])
