#!/usr/bin/env Rscript

library(optparse)
library(Matrix)
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(SingleCellExperiment))

#' Make SingleCellExperiment object from counts, rowData, and colData tables
#'
#' Indexes data by gene common names and spot names.
#' Other processing (e.g. cleaning of row/col columns) should be done upstream
#' in per-dataset scripts
#'
#' @param counts Counts matrix (genes x spots)
#' @param rowData Row data (must include "gene_name" column)
#' @param colData Column data (must include "spot" column)
#' @return SingleCellExperiment object
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

    if ("gene_id" %in% colnames(rowData)) {
        rownames(rowData) <- uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
    } else {
        # TODO: decide how to handle gene names with multiple corresponding IDs
        rownames(rowData) <- rowData$gene_name
    }

    # Make count matrix sparse
    counts <- as(counts, "dgCMatrix")

    sce <- SingleCellExperiment(assays=list(counts=counts),
                                rowData=rowData,
                                colData=colData)
}

add_metadata <- function(sce, options) {
    if ("sample" %in% names(options)) {
        colData(sce)$sample <- options$sample
    }

    return(sce)
}

#' Simple processing of spatial dataset
#'
#' Log-normalized counts + denoised PCA
process_SCE <- function(sce, options) {
    # TODO: benchmark/parameterize the choices here. 
    # Use existing defaults until then

    # Log-normalized counts
    sce <- scater::logNormCounts(sce)

    # De-noised PCA
    dec <- scran::modelGeneVarByPoisson(sce)
    top <- scran::getTopHVGs(dec, prop=0.1)
    sce <- scran::denoisePCA(sce, technical=dec, subset.row=top)

    # Other dimensionality reduction
    # TODO: solve error 
    #   Error in density.default(x, adjust = adjust, from = min(x), to = max(x)) :
    #   need at least 2 points to select a bandwidth automatically
    # sce <- scater::runTSNE(sce, dimred="PCA")
    # sce <- scater::runUMAP(sce, dimred="PCA")
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
    sce <- add_metadata(sce, options)
    sce <- process_SCE(sce, options)

    saveRDS(sce, args[[4]])
}

if (sys.nframe() == 0) {
    main()
}
