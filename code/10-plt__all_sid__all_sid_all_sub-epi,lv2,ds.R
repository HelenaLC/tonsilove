# args <- list(
#     list.files("outs", "^epi", full.names=TRUE),
#     list.files("outs", "lv2", full.names=TRUE),
#     "plts/epi,lv2,ds.pdf")

# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidytext)
    library(HDF5Array)
    library(SingleCellExperiment)
})

.sub <- \(.) gsub(".*([a-z]{3})\\..*", "\\1", .)
.sid <- \(.) split(., gsub(".*([A-Z][0-9]).*", "\\1", .))
args[-3] <- lapply(args[-3], \(.) grep("C1", ., value=TRUE))

df <- mapply(
    x=.sid(args[[1]]),
    y=.sid(args[[2]]),
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        sce <- readRDS(x)
        sid <- paste(sce$sid[1])
        ist <- lapply(y, readRDS)
        kid <- lapply(ist, \(.) .$clust)
        ns <- sapply(kid, length)
        idx <- names(kid <- unlist(kid))
        sub <- rep.int(.sub(y), ns)
        sub <- setNames(sub, idx)
        sub <- sub[match(colnames(sce), idx)]
        kid <- kid[match(colnames(sce), idx)]
        data.frame(sub, kid, colData(sce))
    }) |> do.call(what=rbind)

# wrangling
fd <- df |>
    filter(!is.na(roi)) |>
    mutate(foo=paste(sid, roi, sep=".")) |>
    filter(foo == (id <- "C1.EPI03")) 

# for each subpopulation, 
# estimate distance to apical
mm <- grep("global_mm$", names(fd))
ap <- grepl("epi.api", fd$kid)
xy <- fd[ap, mm]; yx <- fd[!ap, mm]
nn <- nn2(xy, yx, k=10)
ds <- rowMeans(nn$nn.dist)
fd <- data.frame(fd[!ap, ], ds)
fd <- filter(fd, !grepl("epi", kid))

# restrict to most frequent subpopulations
ns <- tail(sort(table(fd$kid)), 20)
fd <- filter(fd, kid %in% names(ns))

# plotting
gg <- ggplot(fd, aes(ds, reorder(kid, ds, median), fill=sub)) +
    scale_fill_manual(NULL, drop=FALSE, values=.pal_sub) +
    labs(x="distance to apical (mm)") +
    ggtitle(.lab(id, nrow(fd))) +
    geom_boxplot(
        alpha=2/3, key_glyph="point", outlier.stroke=0.1,
        outlier.shape=16, outlier.size=0.1, linewidth=0.1) + 
    scale_y_discrete(position="right", limits=\(.) rev(.)) +
    .thm_fig_d("bw", "f") + theme(
        aspect.ratio=2/3,
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.x=element_blank()) 

# saving
ggsave(args[[3]], gg, units="cm", width=5, height=3)
