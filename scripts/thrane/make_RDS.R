#!/usr/bin/env Rscript

library(optparse)
library(Matrix)
suppressMessages(library(scater))
suppressMessages(library(SingleCellExperiment))


make_SCE <- function(counts, rowData, colData) {
    # Load cleaned CSVs
    rowData <- read.csv(rowData)
    colData <- read.csv(colData)
    counts <- as.matrix(read.csv(counts, row.names=1, check.names=F))

    # Index by spot and gene name
	rownames(colData) <- colData$spot
	rownames(rowData) <- uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)

    # Make count matrix sparse
	counts <- as(counts, "dgCMatrix")

    sce <- SingleCellExperiment(assays=list(counts=counts),
                                rowData=rowData,
                                colData=colData)
}

main <- function() {
    # TODO clean this parsing up so it's not as hard-coded
    option_list <- list(
        make_option("--sample")
    )
    parser <- OptionParser(usage = "%prog [options] counts rowData colData RDS",
                           option_list=option_list)
    opt <- parse_args(parser, positional_arguments=4)
    args <- opt$args
    options <- opt$options

    sce <- make_SCE(args[[1]], args[[2]], args[[3]])

    # Add sample label if provided
    if (!is.null(options$sample)) {
        colData(sce)$sample <- options$sample
    }

    saveRDS(sce, args[[4]])
}

if (sys.nframe() == 0) {
    main()
}
