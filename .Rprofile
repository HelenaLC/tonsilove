# dependencies
suppressPackageStartupMessages({
    library(data.table)
    library(BiocParallel)
}); source("code/_utils.R")

# unpack wildcards as named list
.get_wcs <- function(wcs) {
    wcs <- gsub("(,)(\\w+=)", ";\\2", wcs)
    ss <- strsplit(wcs, ";")[[1]]
    ss <- sapply(ss, strsplit, "=")
    keys <- sapply(ss, .subset, 1)
    vals <- sapply(ss, .subset, 2)
	wcs <- as.list(vals)
	names(wcs) <- keys
    return(wcs)
}

# get command line arguments
args <- R.utils::commandArgs(
    trailingOnly = TRUE, 
    asValues = TRUE)

# wildcards
wcs <- args$wcs
if (!is.null(wcs)) {
    wcs <- .get_wcs(wcs)
    args$wcs <- NULL
} else {
    wcs <- NULL
}

# threads
ths <- args$ths
if (!is.null(ths)) {
    args$ths <- NULL
    ths <- as.integer(ths)
    bp <- MulticoreParam(ths)
    setDTthreads(ths)
} else {
    bp <- SerialParam()
    setDTthreads(1)
}

# unpack collapsed arguments
args <- lapply(args, function(u) unlist(strsplit(u, ";")))

# print wildcards & arguments to console
cat("WILDCARDS:\n\n"); print(wcs); cat("\n")
cat("ARGUMENTS:\n\n"); print(args); cat("\n")
