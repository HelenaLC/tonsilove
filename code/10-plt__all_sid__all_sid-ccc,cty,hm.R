# args <- list(
#     list.files("outs", "ccc", full.names=TRUE),
#     list.files("outs", "cty", full.names=TRUE),
#     "plts/ccc,cty,hm.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scuttle)
    library(HDF5Array)
    library(SingleCellExperiment)
})

se <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        ccc <- readRDS(x)
        ctx <- readRDS(y)
        # construct 'SummarizedExperiment'
        cs <- intersect(rownames(ccc[[1]]), colnames(ctx))
        cd <- colData(ctx)[match(cs, colnames(ctx)), ]; rownames(cd) <- cs
        es <- lapply(ccc, \(.)
            `rownames<-`(mx <- t(.[cs, ]),
            gsub("^(s|r)-", "", rownames(mx))))
        # characterize sender/receiver as
        # total, ligand-receptor, pathway
        rd <- data.frame(
            row.names=nm <- rownames(es[[1]]),
            typ=ifelse(grepl("total", nm), "tl",
                ifelse(grepl("-", nm), "lr", "pw")))
        se <- SummarizedExperiment(es, rowData=rd, colData=cd)
    }) |> do.call(what=cbind)

# average by context
names(sr) <- sr <- c("s", "r")
sf <- sg <- se; sf$sid <- ""; sf$foo <- sf$ctx
sg$foo <- paste(sg$ctx, se$sid, sep="-")
sg <- cbind(sf, sg); table(sg$foo)
mu <- lapply(sr, \(.) {
    id <- colData(sg)["foo"]
    sf <- summarizeAssayByGroup(sg, id, assay.type=., statistics="mean")
    df <- data.frame(t(assay(sf)), colData(sf), check.names=FALSE)
    fd <- pivot_longer(df, all_of(rownames(sf)), names_to="lr")
    fd |> mutate(typ=rowData(sg)[fd$lr, "typ"])
}) |> bind_rows(.id="sr") |>
    # average sender/receiver
    group_by(foo, typ, lr) |>
    summarise_at("value", mean) |>
    ungroup()

# get identifiers
mu <- mutate(mu, 
    ctx=gsub("-.*", "", foo),
    sid=gsub(".*-", "", foo),
    sid=case_when(sid %in% names(.pal_sid) ~ sid, TRUE ~ "all"))

# get inter-niche FCs
is <- split(seq(nrow(mu)), mu$sid)
mu <- lapply(is, \(.) {
    mv <- select(mu[., ], c(lr, ctx, value))
    mx <- pivot_wider(mv, names_from="lr")
    my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
    fc <- sapply(rownames(my), \(.) {
        mz <- my[-match(., rownames(my)), ]
        my[., ]/colMeans(mz)
    })
    fc <- data.frame(lr=rownames(fc), fc, row.names=NULL)
    fc <- pivot_longer(fc, -lr, names_to="ctx", values_to="fc")
    left_join(mu[., ], fc, by=c("lr", "ctx"))
}) |> do.call(what=rbind)
mu <- mutate(mu, 
    ctx=factor(ctx, unique(ctx)),
    sid=factor(sid, c("all", names(.pal_sid))))

# plotting
cs <- levels(mu$ctx)
ys <- c("pw", "lr", cs)
ps <- by(mu, mu$sid, \(mv) {
    lapply(ys, \(id) {
        if (id == "pw") {
            n <- 60 
            mv <- mv |> 
                filter(typ == "pw") |> 
                select(-typ)
        } else if (id == "lr") { 
            n <- 15
            mv <- mv |> 
                filter(typ == "lr") |> 
                select(-typ)
        } else { 
            n <- 30 
            cs <- id
            mv <- mv |> 
                filter(typ == "lr") |> 
                select(-typ)
        }
        # selection
        nm <- mv |>
            filter(ctx %in% cs) |>
            group_by(ctx) |>
            slice_max(fc, n=n, with_ties=FALSE) |>
            pull(lr)
        mv <- mv |>
            select(-fc) |>
            filter(lr %in% nm) |>
            group_by(lr) |>
            mutate_at("value", .z)
        # hierarchical clustering
        mx <- pivot_wider(select(mv, c(lr, ctx, value)), names_from="lr")
        my <- `rownames<-`(as.matrix(mx[, -1]), mx[[1]])
        # plotting
        id <- paste(id, mv$sid[1], sep="-")
        ggplot(mv, aes(ctx, lr, fill=value)) + 
            geom_tile() + ggtitle(id) +
            scale_fill_gradient2(
                "z-scaled\nmean CCC",
                limits=c(-2.5, 2.5), n.breaks=6,
                low="cadetblue", mid="ivory", high="firebrick") +
            scale_x_discrete(limits=\(.) intersect(.xo(my), .)) +
            scale_y_discrete(limits=.yo(my)) +
            coord_equal(2/3, expand=FALSE) +
            .thm_fig_c("minimal") + theme(
                panel.grid=element_blank(),
                axis.title=element_blank(),
                axis.text.y=element_text(size=2),
                axis.text.x=element_text(size=3, angle=90, hjust=1, vjust=0.5)) 
    }) 
}) |> Reduce(f=c)

# saving
tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
for (. in seq_along(ps)) {
    x <- ps[[.]]$data
    w <- nlevels(x$ctx)/8+max(nchar(x$lr))/12
    h <- length(unique(x$lr))/12+max(nchar(levels(x$ctx)))/12
    pdf(tf[[.]], 
        width=(1+w)/2.54, 
        height=(0.5+h)/2.54)
    print(ps[[.]]); dev.off()
}
qpdf::pdf_combine(unlist(tf), output=args[[3]])
