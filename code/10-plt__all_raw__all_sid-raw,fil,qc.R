# args <- list(
#     list.files("outs", "se\\.rds", recursive=TRUE, full.names=TRUE),
#     list.files("outs", "fil", full.names=TRUE),
#     "plts/raw,fil,qc.pdf")

# dependencies
suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# loading
sce <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        raw <- readRDS(x)
        fil <- readRDS(y)
        sid <- as.character(raw$sid[1])
        assays(raw) <- altExps(raw) <- list()
        raw$fil <- colnames(raw) %in% colnames(fil)
        metadata(raw)[[sid]] <- metadata(fil)$ths
        raw
    }) |> do.call(what=cbind)

# wrangling
xs <- c("counts", "features", "area", "counts/um2", "features/um2")
cd <- data.frame(colData(sce))
df <- cd |>
    mutate(area=Area.um2, 
        counts=nCount_RNA, 
        features=nFeature_RNA) |>
    mutate(`counts/um2`=counts/area) |>
    mutate(`features/um2`=features/area) |>
    mutate(sid=factor(sid, unique(sid))) |>
    mutate(fil=factor(fil,
        levels=c(TRUE, FALSE),
        labels=c("kept", "removed"))) |>
    pivot_longer(all_of(xs)) |>
    mutate(name=factor(name, xs)) |>
    mutate(sid=factor(sid, c("C1", "A1", "A2")))

# annotations
mu <- df |>
    group_by(sid, name) |>
    summarise_at("value", mean) |>
    mutate(lab=round(value, 2)) |>
    mutate(lab=format(lab, big.mark=","))

th <- do.call(cbind, lapply(metadata(sce), sapply, as.numeric))
th <- data.frame(name=rownames(th), th)
nm <- c("counts", "features", "counts/um2")
th$name <- nm[match(th$name, c("n", "m", "a"))]
th$name <- factor(th$name, levels(df$name))
th <- pivot_longer(th, -name, names_to="sid") 
th <- th |>
    filter(!is.na(name)) |>
    mutate(lab=round(value, 2)) |>
    mutate(lab=format(lab, big.mark=",")) |>
    mutate(sid=factor(sid, levels(df$sid))) |>
    mutate(name=factor(name, levels(df$name))) 

# plotting
gg <- ggplot(df, aes(value, after_stat(ndensity))) + facet_grid(sid~name, scales="free") + 
    geom_histogram(bins=30, fill="blue", col="lavender", alpha=1/3, linewidth=0.1) +
    geom_vline(data=mu, aes(xintercept=value), col="blue", linewidth=0.2) +
    geom_vline(data=th, aes(xintercept=value), col="red", linewidth=0.2) +
    geom_text(data=mu, aes(Inf, Inf, label=lab), vjust=1.2, hjust=1.1, size=1.2, col="blue") +
    geom_text(data=th, aes(Inf, Inf, label=lab), vjust=2.3, hjust=1.1, size=1.2, col="red") +
    scale_x_continuous(NULL, transform="log10") +
    scale_y_continuous("normalized density", n.breaks=3) + 
    .thm_fig("bw") + theme(aspect.ratio=2/3)

# saving
ggsave(args[[3]], gg, unit="cm", width=12, height=6)
