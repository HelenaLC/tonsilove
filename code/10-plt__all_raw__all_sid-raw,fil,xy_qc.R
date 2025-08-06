# args <- list(
#     c("outs/raw-A1", "outs/raw-C1"), 
#     c("outs/fil-A1.rds", "outs/fil-C1.rds"), 
#     "plts/raw,fil,xy_qcs.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})

ps <- lapply(args[[1]], \(.) {
    # loading
    se <- readRDS(.)
    cd <- colData(se)
    df <- data.frame(cd)
    # wrangling
    fd <- df |>
        mutate(area=Area.um2, 
            counts=nCount_RNA, 
            features=nFeature_RNA) |>
        mutate(`counts/um2`=counts/area) |>
        mutate(`features/um2`=features/area)
    # plotting
    xs <- setdiff(names(fd), names(df))
    fd <- mutate(fd, across(setdiff(xs, "area"), log10))
    lapply(xs, \(.) {
        .plt_xy(fd, .z(fd[[.]])) + ggtitle(.) +
            scale_color_gradientn(
                limits=c(-2.5, 2.5),
                breaks=seq(-2, 2, 2),
                colors=hcl.colors(9, "Plasma"))
    })
}) |> Reduce(f=c)

# saving
.pdf(ps, args[[3]])
