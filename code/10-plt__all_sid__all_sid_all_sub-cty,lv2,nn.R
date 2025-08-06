# # args <- list(
# #     list.files("outs", "cty", full.names=TRUE),
# #     list.files("outs", "lv2", full.names=TRUE),
# #     "plts/cty,lv2,nn.pdf")

# dependencies
suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
    library(tidyr)
    library(HDF5Array)
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
df <- by(df, df[c("sid", "ctx")], \(fd) {
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
    cbind(fd[c("sid", "kid", "ctx")], fq)
}) |> do.call(what=rbind)

# plotting
ps <- by(df, df[c("sid", "ctx")], \(fd) {
    sid <- paste(fd$sid[1])
    ctx <- paste(fd$ctx[1])
    id <- paste(ctx, sid, sep="-")
    ks <- tail(names(sort(table(fd$kid))), 15)
    # averaging
    mu <- fd |>
        group_by(kid) |>
        filter(!is.na(kid)) |>
        summarise(
            .groups="drop",
            across(where(is.numeric),
            \(.) mean(., na.rm=TRUE)))
    # wrangling
    mv <- mu |>
        filter(kid %in% ks) |>
        pivot_longer(where(is.numeric)) |>
        filter(name %in% ks) |>
        mutate(value=case_when(kid==name~NA, TRUE~value)) |>
        group_by(name) |> mutate_at("value", .z)
    # hierarchical clustering
    mx <- pivot_wider(mv)
    my <- as.matrix(mx[, -1])
    rownames(my) <- mx[[1]]
    # plotting
    nm <- "z-scaled mean frequency of X\nin radial neighborhood of Y"
    gg <- ggplot(mv, aes(name, kid, fill=value)) + 
        geom_tile() + coord_equal(expand=FALSE) +
        scale_fill_gradient2(nm,
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2),
            low="gold", high="navy", na.value="grey") +
        scale_x_discrete(limits=.xo(my)) +
        scale_y_discrete(limits=.xo(my)) +
        ggtitle(.lab(id, nrow(fd))) +
        .thm_fig_c("minimal") + theme(
            legend.position="bottom",
            legend.title.position="top",
            axis.title=element_blank(),
            legend.key.width=unit(0.8, "lines"),
            legend.key.height=unit(0.2, "lines"),
            legend.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
})
    
# saving
pdf(args[[3]], onefile=TRUE, width=3.5/2.54, height=4/2.54)
for (p in ps) print(p); dev.off()
