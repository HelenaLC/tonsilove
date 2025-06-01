#args <- list(list.files("outs", "ist", full.names=TRUE), "plts/ist,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
names(ist) <- sid <- gsub(".*((A|C)[1-2]).*", "\\1", args[[1]])
df <- do.call(rbind, lapply(sid, \(sid) data.frame(sid, kid=ist[[sid]]$clust)))

# plotting
s <- length(unique(df$sid))
k <- length(unique(df$kid))
p <- .plt_fq(df, "sid", "kid", hc=TRUE, a=0.8) +
    scale_fill_manual("cluster", values=.pal) +
    labs(x=NULL) + theme(aspect.ratio=8) 
q <- .plt_fq(df, "kid", "sid", hc=TRUE) + 
    scale_fill_manual("section", values=.pal_sid) +
    labs(x=NULL) + theme(aspect.ratio=8*s/k)
gg <- wrap_plots(p, q, nrow=1) +
    plot_layout(guides="collect") &
    theme(plot.title=element_blank(),
        legend.spacing=unit(0, "lines")) &
    guides(fill=guide_legend(override.aes=list(
        alpha=1, shape=21, stroke=0, size=2)))

# saving
ggsave(args[[2]], gg, units="cm", width=9, height=6)

