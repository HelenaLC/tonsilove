# args <- list(
#     list.files("outs", "se\\.rds", recursive=TRUE, full.names=TRUE),
#     list.files("outs", "fil", full.names=TRUE),
#     "plts/raw,fil,qc.pdf")

# dependencies
suppressPackageStartupMessages({
    library(zoo)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(SingleCellExperiment)
})

ps <- mapply(
    x=args[[1]], y=args[[2]],
    SIMPLIFY=FALSE, \(x, y) {
        # loading
        raw <- readRDS(x)
        fil <- readRDS(y)
        
        # wrangling
        cd <- data.frame(colData(raw))
        cd$cpa <- with(cd, nCount_RNA/Area.um2)
        df <- cbind(cd, d <- .d2b(cd)) |>
            pivot_longer(all_of(colnames(d))) |>
            mutate(value=value*1e3) |>
            mutate(name=factor(name, 
                c("t", "r", "b", "l"),
                c("top", "right", "bottom", "left"))) |>
            group_by(name, value) |>
            summarise_at("cpa", mean) |>
            arrange(value) |> group_by(name) |> mutate(
                mu_x=rollmean(value, 20, mean, align="left", fill=0),
                mu_y=rollmean(cpa, 20, mean, align="left", fill=0)) 
        mu <- data.frame(name=unique(df["name"]), y=mean(cd$cpa))
        
        # aesthetics
        dy <- ceiling(max(df$mu_y)/0.1)*0.1
        n <- ncol(raw); m <- ncol(fil)
        p <- round(100*m/n, 2)
        n <- format(n, big.mark=",")
        m <- format(m, big.mark=",")
        i <- paste(raw$sid[1])
        md <- metadata(fil)$ths
        
        # plotting
        ggplot(df, aes(mu_x, mu_y)) + 
            facet_grid(~name) + geom_line(col="navy", linewidth=0.2) +
            geom_vline(xintercept=md$d, col="red", linewidth=0.2, lty=2) +
            geom_hline(yintercept=md$a, col="red", linewidth=0.2, lty=2) +
            geom_hline(data=mu, aes(yintercept=y), col="blue", linewidth=0.2) +
            scale_x_continuous("FOV border distance (um)", expand=expansion(0.04)) +
            scale_y_continuous("RNA counts per um", n.breaks=3, expand=expansion(0.06)) +
            ggtitle(bquote(bold(.(i))~"(N ="~.(n)*";"~.(m)*";"~.(p)*"%)")) +
            coord_cartesian(xlim=c(0, 60), ylim=c(0, dy)) +
            .thm_fig("bw") + theme(aspect.ratio=2/3) 
    })

# saving
gs <- gridExtra::marrangeGrob(grobs=ps, nrow=1, ncol=1, top=FALSE)
ggsave(args[[3]], gs, unit="cm", width=10, height=3)

