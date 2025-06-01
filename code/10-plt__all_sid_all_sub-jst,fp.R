# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(InSituType)
})

# loading
ist <- lapply(args[[1]], readRDS)

# plotting
sub <- gsub(".*,([a-z]{3}).*", "\\1", args[[1]])
sid <- gsub(".*((A|C)[0-9]).*", "\\1", args[[1]])
ps <- lapply(split(seq_along(ist), sub), \(i) {
    ps <- lapply(seq_along(ist[i]), \(j) {
        .plt_fp(ist[i][[j]], sid[i][j]) })
    wrap_plots(ps, nrow=1) & 
    theme(plot.margin=margin(l=1, r=1))
})

# saving
pdf(args[[2]], onefile=TRUE, width=12/2.54, height=4.5/2.54)
for (p in ps) print(p); dev.off()
