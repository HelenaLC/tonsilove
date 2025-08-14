foo <- lapply(args[[1]], \(x) {
    # loading
    sce <- readRDS(x)
    # labeling
    lb <- list(
        `B/T`="N5", TCZ="N2", 
        CTS=c("N7", "N10"), 
        EPI="N3", SM="N4", 
        GC=c("N1", "N8", "N9"), MZ="N6")
    cs <- match(sce$ctx, unlist(lb))
    lb <- rep.int(names(lb), sapply(lb, length))
    table(sce$ctx, sce$ctx <- lb[cs]); table(sce$ctx)
    # saving
    y <- sub("ctx", "cty", x)
    saveRDS(sce, y)
})
