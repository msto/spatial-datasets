#!/usr/bin/env Rscript

library(optparse)
library(Matrix)
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(SingleCellExperiment))

make_SCE <- function(counts, rowData, colData) {
    # Load cleaned CSVs
    rowData <- read.csv(rowData)
    colData <- read.csv(colData)
    counts <- as.matrix(read.csv(counts, row.names=1, check.names=F))

    if (dim(counts)[1] != dim(rowData)[1]) {
        stop("Count matrix and rowData contain different number of genes")
    }

    if (dim(counts)[2] != dim(colData)[1]) {
        stop("Count matrix and colData contain different number of spots")
    }

    if (!("spot" %in% colnames(colData))) {
        stop("colData is missing 'spot' ID column")
    }

    if (!("gene_name" %in% colnames(rowData))) {
        stop("rowData is missing 'gene_name' ID column")
    }

    # TODO: check rowname ordering across counts and row/coldata

    # Index by spot and gene name
    rownames(colData) <- colData$spot

    # TODO: decide how to handle gene names with multiple corresponding IDs
    # rownames(rowData) <- uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
    rownames(rowData) <- rowData$gene_name

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

    sce <- scater::logNormCounts(sce)
    # sce <- scran::denoisePCA(sce)

    saveRDS(sce, args[[4]])
}

if (sys.nframe() == 0) {
    main()
}
