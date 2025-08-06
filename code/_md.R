sce <- readRDS("outs/fil-A2.rds")
ctx <- readRDS("outs/cty-A2.rds")
ist <- readRDS("outs/lv1-A2.rds")$clust
sce$ist <- ist[match(colnames(sce), names(ist))]
sce$ctx <- ctx$ctx[match(colnames(sce), colnames(ctx))]
df <- data.frame(cell_ID=sce$cell, colData(sce))
write.table(df, "imgs/S2/metadata.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
