#args <- list(list.files("outs", "lv1", full.names=TRUE), "plts/lv1,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
names(ist) <- sid <- gsub(".*((A|C)[1-2]).*", "\\1", args[[1]])
df <- do.call(rbind, lapply(sid, \(sid) data.frame(sid, kid=ist[[sid]]$clust)))
ns <- arrange(dplyr::count(group_by(df, kid)), n) |> mutate(m=format(n, big.mark=","))
ns$n[ns$n > 2e5] <- 2e5

# plotting
k <- length(unique(df$kid))
s <- length(unique(df$sid))
n <- ggplot(ns, aes(n, kid, fill=kid)) + 
    geom_bar(stat="identity", show.legend=FALSE) +
    geom_text(aes(x=2e5, label=m), size=1.5, hjust=1) +
    scale_fill_manual(values=.pal) +
    scale_x_continuous(
        breaks=c(0, 1e5, 2e5),
        labels=c(0, "100k", ">200k")) +
    coord_cartesian(expand=FALSE) +
    .thm_fig("classic") + theme(
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        aspect.ratio=1)

p <- .plt_fq(df, "sid", "kid", h=TRUE) +
    scale_fill_manual("cluster", values=.pal) +
    scale_x_discrete(limits=c("A2", "A1", "C1")) +
    theme(aspect.ratio=s/k)

q <- .plt_fq(df, "kid", "sid", h=TRUE) + 
    scale_fill_manual("section", values=.pal_sid) +
    scale_x_discrete(limits=ns$kid) +
    theme(aspect.ratio=1)

o <- q$scales$scales[[2]]$limits
n <- n + scale_y_discrete(limits=o)

gg <- wrap_plots(nrow=1,
    wrap_plots(p, q, ncol=1, heights=hs <- c(s, k)),
    wrap_plots(plot_spacer(), n, ncol=1, heights=hs)) +
    plot_layout(guides="collect") &
    theme(
        plot.title=element_blank(),
        axis.title=element_blank(),
        legend.spacing=unit(0, "lines")) &
    guides(fill=guide_legend(override.aes=list(
        alpha=1, shape=21, stroke=0, size=1)))

# saving
ggsave(args[[2]], gg, units="cm", width=9, height=5)

