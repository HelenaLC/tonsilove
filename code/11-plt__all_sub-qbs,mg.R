# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(scater)
    library(scuttle)
    library(ggplot2)
    library(SingleCellExperiment)
})

# loading
sub <- gsub(".*-([a-z]{3}).*", "\\1", args[[1]])
lys <- lapply(setNames(args[[1]], sub), readRDS)

# automated selection
mgs <- lapply(lys, \(sce) {
    sce <- sce[, sce$sid == "all"]
    colnames(sce) <- names(ks) <- ks <- sce$kid
    es <- logcounts(sce)
    gs <- lapply(ks, \(i) {
        j <- setdiff(ks, i)
        fc <- es[, i]/rowMeans(es[, j])
        names(tail(sort(fc), 10))
    })
    is <- split(seq_along(unlist(gs)), rep(names(gs), sapply(gs, length)))
    is <- lapply(is, \(.) setdiff(., which(duplicated(unlist(gs)))))
    lapply(is, \(.) unlist(gs)[.])
})

# curated selection
# mgs <- list(
#     bcs=list(
#         
#     ),
#     tcs=list(
#         
#     ),
#     mye=list(
#         
#     ),
#     str=list(
#         
#     )
# )
lapply(mgs, \(.) unlist(.)[duplicated(unlist(.))])

# aesthetics
aes <- list(
    geom_tile(), 
    coord_equal(4/3, expand=FALSE),
    .thm_fig_c("minimal"), theme(
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5)),
    scale_fill_gradient2(
        "z-scaled\nmean expr.", 
        low="dodgerblue2", high="tomato2",
        limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)))

ps <- lapply(sub, \(.) {
    se <- lys[[.]]
    gs <- mgs[[.]]
    ns <- cumsum(sapply(gs, length))
    ns <- unique(ns); ns <- ns[-length(ns)]
    # joint
    x <- se[, se$sid == "all"]
    y <- t(logcounts(x[unlist(gs), ]))
    df <- data.frame(y, k=x$kid, check.names=FALSE)
    fd <- df |>
        pivot_longer(-k, names_to="g", values_to="y") |>
        group_by(g) |> mutate_at("y", .z) |>
        mutate_at("k", factor, rev(names(gs))) |>
        mutate_at("g", factor, unlist(gs))
    p <- ggplot(fd, aes(g, k, fill=y)) + aes + ggtitle(.) +
        geom_vline(xintercept=ns+0.5, linewidth=0.2)
    # split
    x <- se[, se$sid != "all"]
    y <- t(logcounts(x[unlist(gs), ]))
    df <- data.frame(y, s=x$sid, k=x$kid, check.names=FALSE)
    df <- mutate(df, across(where(is.numeric), .z))
    q <- by(df, df$k, \(fd) {
        k <- fd$k[1]
        fd <- fd |>
            select(-k) |>
            pivot_longer(-s, names_to="g", values_to="y") |>
            mutate_at("g", factor, unlist(gs))
        # plotting
        hl <- c("plain", "bold")[1 + fd$g %in% gs[[k]]]
        ggplot(fd, aes(g, s, fill=y)) + aes + ggtitle(k) +
            geom_vline(xintercept=ns+0.5, linewidth=0.2) +
            scale_y_discrete(limits=\(.) rev(.)) + guides(fill="none") +
            theme(axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5, face=hl))
    })
    list(p)#c(list(p), q)
}) |> Reduce(f=c)

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    df <- (p <- ps[[.]])$data
    y <- ifelse(is.null(df$k), "s", "k")
    w <- length(unique(df$g))/10+max(nchar(paste(df[[y]])))/5
    h <- length(unique(df[[y]]))/10+max(nchar(paste(df$g)))/5
    pdf(tf[[.]], width=(2+w)/2.54, height=(0.5+h)/2.54)
    print(p); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[2]])
