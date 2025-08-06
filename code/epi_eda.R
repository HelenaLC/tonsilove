args <- list(
    list.files("outs", "^epi", full.names=TRUE),
    list.files("outs", "lv2", full.names=TRUE),
    "plts/epi,lv2,pca.pdf")

# # dependencies
# suppressPackageStartupMessages({
#     library(dplyr)
#     library(ggrepel)
#     library(HDF5Array)
#     library(SpatialExperiment)
#     library(SingleCellExperiment)
# })
# 
# .sub <- \(.) gsub(".*([a-z]{3})\\..*", "\\1", .)
# .sid <- \(.) split(., gsub(".*([A-Z][0-9]).*", "\\1", .))
# se <- mapply(
#     x=.sid(args[[1]]),
#     y=.sid(args[[2]]),
#     SIMPLIFY=FALSE, \(x, y) {
#         # loading
#         sce <- readRDS(x)
#         sid <- paste(sce$sid[1])
#         ist <- lapply(y, readRDS)
#         # wrangling
#         mtx <- ist[[grep("epi", y)]]$profiles
#         sel <- rownames(sce) %in% rownames(mtx)
#         rowData(sce)$sel <- sel
#         kid <- lapply(ist, \(.) .$clust)
#         ns <- sapply(kid, length)
#         idx <- names(kid <- unlist(kid))
#         sub <- rep.int(.sub(y), ns)
#         sub <- setNames(sub, idx)
#         sub <- sub[match(colnames(sce), idx)]
#         kid <- kid[match(colnames(sce), idx)]
#         cd <- cbind(colData(sce), data.frame(sub, kid))
#         `colData<-`(sce, value=cd)
#     })
# se <- lapply(se, as, "SingleCellExperiment")
# (se <- do.call(cbind, se))
# 
# # aggregation
# sf <- se[, !is.na(se$roi) & se$sub == "epi"]
# sf$roj <- paste(sf$sid, sf$roi, sep=".")
# id <- colData(sf)[c("roj")]
# sf <- logNormCounts(sf)
gs <- rowData(sf)$sel
# 
# pb <- aggregateAcrossCells(sf, id,
#     use.assay.type="logcounts",
#     statistics="mean")

# specify substructures that correspond to surface epithelium
sq <- c(
    "A1.EPI01", "A1.EPI05", "A2.EPI19", "A2.EPI23",
    "A2.EPI03", "A2.EPI04", "A2.EPI07", "A2.EPI08", 
    "A2.EPI09", "A2.EPI12", "A2.EPI13", "A2.EPI20", 
    "A2.EPI24", "A2.EPI26", "A2.EPI27")
pb$typ <- ifelse(pb$roj %in% sq, "sq", "cr")

gs <- grep("KRT", rownames(pb), value=TRUE)
pal <- pals::polychrome()[-2]
pal <- setNames(pal[seq_along(rs)], rs <- levels(pb$roi))
plotHeatmap(pb[, pb$sid == "A2"], gs, 
    center=TRUE,
    color_columns_by=c("sid", "roi", "typ"),
    column_annotation_colors=list(
        sid=.pal_sid, roi=pal, 
        typ=c(sq="black", cr="grey")))

ps <- lapply(c("A1", "A2"), \(.) {
    qb <- pb[, pb$sid == .]
qb <- scater::runPCA(qb, subset_row=gs, ncomponents=5)

# plotting
pc <- reducedDim(qb, "PCA")
df <- data.frame(colData(qb), pc)
df$sid <- factor(df$sid, names(.pal_sid))

rot <- attr(reducedDim(qb, "PCA"), "rotation")
top <- names(tail(sort(rowMaxs(abs(rot))), 40))
rot <- data.frame(var=rownames(rot), rot[, 1:2])[top, ]

# plot(seq(ncol(pc)), attr(pc, "varExplained"), type="b")
# ks <- kmeans(pc[, seq(6)], centers=3)$clust
# 
# ggplot(mutate(df, k=factor(ks)), aes(PC1, PC2, fill=k, shape=typ)) + 
#     geom_vline(xintercept=0, linewidth=0.1) +
#     geom_hline(yintercept=0, linewidth=0.1) +
#     scale_shape_manual(values=c(21, 24)) +
#     #scale_fill_manual(values=.pal_sid) +
#     geom_point(size=2, stroke=0) +
#     geom_text_repel(
#         aes(label=gsub("EPI", "", roi)), size=1.2,
#         min.segment.length=0, segment.size=0.1) +
#     .thm_fig_d("bw", "f") + theme()

# components
pal <- setNames(pals::polychrome()[-2][seq_along(rs)], rs <- levels(df$roi))
p <- ggplot(df, aes(PC1, PC2, fill=roi, shape=typ)) +
    geom_vline(xintercept=0, linewidth=0.1) +
    geom_hline(yintercept=0, linewidth=0.1) +
    scale_shape_manual(values=c(21, 24)) +
    scale_fill_manual(values=pal) +
    geom_point(size=2, stroke=0) +
    geom_text_repel(
        aes(label=gsub("EPI", "", roi)), size=1,
        min.segment.length=0, segment.size=0.1) +
    .thm_fig_d("bw", "f") 
p$guides$guides$fill$params$override.aes$size <- 1

# loadings
q <- ggplot(rot[top, ], aes(PC1, PC2)) +
    geom_vline(xintercept=0, linewidth=0.1) +
    geom_hline(yintercept=0, linewidth=0.1) +
    geom_segment(
        linewidth=0.1, col="grey",
        aes(0, 0, xend=PC1*4, yend=PC2*4),
        arrow=arrow(length=unit(2, "pt"))) +
    annotate("text", size=1, x=rot$PC1*4.5, y=rot$PC2*4.5, label=rot$var) +
    .thm_fig_d("bw", "f") + theme()

gg <- wrap_plots(q, p, nrow=1)
})

# saving
pdf(args[[3]], onefile=TRUE, width=12/2.54, height=5/2.54)
for (p in ps) print(p); dev.off()
