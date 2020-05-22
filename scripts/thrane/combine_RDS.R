#!/usr/bin/env Rscript

library(optparse)
library(Matrix)
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(purrr))


main <- function() {
    option_list <- list()
    parser <- OptionParser(usage = "%prog [options] samples RDS",
                           option_list=option_list)
    opt <- parse_args(parser, positional_arguments=2)
    args <- opt$args
    options <- opt$options

    samples <- scan(args[[1]], character(), quote="")
    sce_vec <- vector("list", length(samples))

    i <- 1
    for (sample in samples) {
        sce_vec[[i]] <- readRDS(sprintf("data/cleaned/thrane/%s.rds", sample))
        i <- i + 1
    }

    genelists <- map(sce_vec, rownames)
    universe <- reduce(genelists, intersect)
    filtered <- map(sce_vec, function(x) x[universe, ])

    # Can't directly cbind SCE objects because of redundant row/colData (e.g. gene names)
    # Until we have a better solution, I'm just cbinding the counts matrices and 
    # then creating a new SCE with updated row/colData
    counts <- reduce(map(filtered, SingleCellExperiment::counts), cbind)
    rowData <- SingleCellExperiment::rowData(filtered[[1]])
    colData <- reduce(map(filtered, SingleCellExperiment::colData), rbind)

    sce <- SingleCellExperiment(assays=list(counts=counts),
                                rowData=rowData, colData=colData)

    saveRDS(sce, args[[2]])
}

if (sys.nframe() == 0) {
    main()
}
