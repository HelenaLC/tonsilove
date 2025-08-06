#args <- list(list.files("outs", "^epi", full.names=TRUE), "plts/epi,xy.pdf")

ps <- lapply(args[[1]], \(x) {
    # loading
    sce <- readRDS(x)
    sid <- paste(sce$sid[1])
    # plotting
    foo <- setNames(sce$roi, colnames(sce))
    .plt_xy(sce, foo, sid, na=TRUE, split=FALSE) + 
        scale_color_manual(
            NULL, na.value="lightgrey", 
            breaks=\(.) setdiff(., NA), 
            values=unname(pals::polychrome()[-2]))
})

# saving
.pdf(ps, args[[2]])
