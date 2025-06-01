# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(HDF5Array)
})

p <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
    # loading
    raw <- readRDS(x)
    fil <- readRDS(y)
    # wrangling
    ex <- !colnames(raw) %in% colnames(fil)
    ex <- setNames(ex, colnames(raw))
    df <- data.frame(colData(raw), ex)
    df <- df[order(df$ex), ]
    # aesthetics
    i <- paste(raw$sid[1])
    n <- ncol(raw)-ncol(fil)
    n <- format(n, big.mark=",")
    p <- round(100*mean(!ex), 2)
    # plotting
    .plt_xy(df, df$ex, split=FALSE) +
        theme(legend.position="none") +
        scale_color_manual(values=c("lavender", "blue")) +
        ggtitle(bquote(bold(.(i))~"(N ="~.(n)*";"~.(p)*"%)"))
})

# saving
.pdf(p, args[[3]])
