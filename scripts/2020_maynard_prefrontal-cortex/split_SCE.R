#!/usr/bin/env Rscript

library(optparse)
suppressMessages({
    library(BayesSpace)
})

main <- function() {
    parser <- OptionParser(usage = "%prog [options] sample RDS",
                           option_list=list())
    opt <- parse_args(parser, positional_arguments=2)
    args <- opt$args

    sample <- args[[1]]
    sce <- spatialLIBD::fetch_data(type="sce")
    sce <- sce[, sce$sample_name == sample]

    sce <- spatialPreprocess(sce, n.PCs=15, n.HVGs=2000)

    saveRDS(sce, args[[2]])
}

if (sys.nframe() == 0) {
    main()
}
