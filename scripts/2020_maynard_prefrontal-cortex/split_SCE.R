#!/usr/bin/env Rscript

library(optparse)
suppressMessages(library(SingleCellExperiment))

main <- function() {
    parser <- OptionParser(usage = "%prog [options] sample RDS",
                           option_list=list())
    opt <- parse_args(parser, positional_arguments=2)
    args <- opt$args

    sample <- args[[1]]
    sce <- spatialLIBD::fetch_data(type="sce")
    sce <- sce[, sce$sample_name == sample]

    set.seed(149)
    dec <- scran::modelGeneVar(sce)
    top <- scran::getTopHVGs(dec, n=2000)
    sce <- scater::runPCA(sce, subset_row=top, ncomponents=15)

    saveRDS(sce, args[[2]])
}

if (sys.nframe() == 0) {
    main()
}
