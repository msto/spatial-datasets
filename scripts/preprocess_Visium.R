#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(BayesSpace)
})


main <- function() {
    args <- commandArgs(trailingOnly=TRUE)

    sce <- readVisium(args[[1]])
    sce <- sce[, colSums(counts(sce)) > 0]

    set.seed(101)
    sce <- spatialPreprocess(sce)

    saveRDS(sce, args[[2]])
}

if (sys.nframe() == 0) {
    main()
}
