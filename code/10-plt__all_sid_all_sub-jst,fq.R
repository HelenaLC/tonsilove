#args <- list(list.files("outs", "jst", full.names=TRUE), "plts/jst,fq.pdf")

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# loading
ist <- lapply(args[[1]], readRDS)

# wrangling
df <- mapply(
    sid=gsub("jst-(.*),.*\\.rds", "\\1", basename(args[[1]])),
    sub=gsub("jst-.*,(.*)\\.rds", "\\1", basename(args[[1]])),
    kid=lapply(ist, \(lys) unname(lys$clust)), SIMPLIFY=FALSE,
    \(sid, sub, kid) data.frame(sid, sub, kid)) |>
    do.call(what=rbind) 

# plotting
is <- split(seq(nrow(df)), df$sub)
gg <- lapply(names(is), \(.) {
    fd <- df[is[[.]], ]
    fd$kid <- as.character(fd$kid)
    .plt_fq(fd, "sid", "kid", id=., hc=FALSE) +
        scale_x_discrete(limits=rev(c("C1", "A1", "A2"))) +
        ggtitle(.lab(., nrow(fd), TRUE))
}) |> wrap_plots(nrow=1) & theme(
    axis.title=element_blank(),
    legend.title=element_blank())

# saving
ggsave(args[[2]], gg, units="cm", width=15, height=6)
