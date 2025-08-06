# args <- list(
#     list.files("outs", "fil", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/fil,lv2,pck.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidytext)
    library(HDF5Array)
    library(SingleCellExperiment)
})

# loading
.f <- \(.) split(., gsub(".*([A-Z][0-9]).*", "\\1", .))
df <- mapply(
    x=.f(args[[1]]),
    y=.f(args[[2]]),
    SIMPLIFY=FALSE, \(x, y) {
        sce <- readRDS(x)
        ist <- lapply(y, readRDS)
        kid <- unlist(lapply(ist, \(.) .$clust))
        kid <- kid[match(colnames(sce), names(kid))]
        data.frame(colData(sce), kid)
    }) |> do.call(what=rbind)

# filtering
fd <- filter(df, grepl("^epi", kid))
fq <- prop.table(table(fd$kid, fd$sid), 2)
ks <- apply(fq, 2, \(.) rownames(fq)[. > 0], simplify=FALSE)
ks <- paste(rep.int(names(ks), sapply(ks, length)), unlist(ks), sep=".")

# wrangling
fd <- fd |>
    mutate(foo=paste(sid, kid, sep=".")) |>
    filter(foo %in% ks) |>
    mutate(pck=asinh(Mean.PanCK/200)) |>
    group_by(sid) |> mutate(pck=.z(pck)) |>
    mutate(sid=factor(sid, names(.pal_sid)))

# plotting
gg <- ggplot(fd, aes(reorder(kid, pck, median), pck, group=foo, fill=sid)) + 
    scale_fill_manual("section", values=.pal_sid) +
    scale_x_reordered(NULL) + labs(y="PanCK") +
    geom_boxplot(
        alpha=2/3, key_glyph="point", outlier.stroke=0.1,
        outlier.shape=16, outlier.size=0.1, linewidth=0.1) + 
    geom_hline(yintercept=0, linewidth=0.1) +
    .thm_fig_d("bw", "f") + theme(
        aspect.ratio=1,
        panel.grid.major.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))

# saving
ggsave(args[[3]], gg, units="cm", width=4, height=2.5)
