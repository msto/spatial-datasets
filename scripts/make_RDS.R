#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(Matrix)
    library(SingleCellExperiment)
    library(BayesSpace)
})

make_SCE_from_10X <- function(dirname) {
    spatial_dir <- file.path(dirname, "spatial")
    matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")

    if (!dir.exists(matrix_dir)) stop(paste0("Matrix directory does not exist: ", matrix_dir))
    if (!dir.exists(spatial_dir)) stop(paste0("Spatial directory does not exist: ", spatial_dir))

    colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header=F)
    # colnames(colData) <- c("spot", "in_tissue", "x", "y", "image_x", "image_y")
    colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
    rownames(colData) <- colData$spot
    colData <- colData[colData$in_tissue > 0, ]

    rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"), header=F)
    colnames(rowData) <- c("gene_id", "gene_name", "X")
    rowData <- rowData[, c("gene_id", "gene_name")]
    rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)

    counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
    barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), header=F)
    colnames(counts) <- barcodes$V1
    rownames(counts) <- rownames(rowData)
    counts <- counts[, rownames(colData)]
    
    sce <- SingleCellExperiment(assays=list(counts=counts),
                                rowData=rowData,
                                colData=colData)
}

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
make_SCE <- function(counts, rowData, colData, options) {
    # Load cleaned CSVs
    rowData <- read.csv(rowData, stringsAsFactors=FALSE)
    colData <- read.csv(colData, stringsAsFactors=FALSE, row.names=1)
    counts <- as.matrix(read.csv(counts, row.names=1, check.names=F, stringsAsFactors=FALSE))

    if (dim(counts)[1] != dim(rowData)[1]) {
        stop("Count matrix and rowData contain different number of genes")
    }

    if (dim(counts)[2] != dim(colData)[1]) {
        stop("Count matrix and colData contain different number of spots")
    }

    if (!("gene_name" %in% colnames(rowData))) {
        stop("rowData is missing 'gene_name' ID column")
    }

    # TODO: check rowname ordering across counts and row/coldata

    if ("gene_id" %in% colnames(rowData)) {
        rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
    } else {
        # TODO: decide how to handle gene names with multiple corresponding IDs
        rownames(rowData) <- rowData$gene_name
    }

    # Make count matrix sparse
    counts <- as(counts, "dgCMatrix")

    if (options$logcounts) {
        sce <- SingleCellExperiment(assays=list(logcounts=counts),
                                    rowData=rowData,
                                    colData=colData)
    } else {
        sce <- SingleCellExperiment(assays=list(counts=counts),
                                    rowData=rowData,
                                    colData=colData)
    }

    sce
}

add_metadata <- function(sce, options) {
    for (key in c("sample", "dataset")) {
        if (key %in% names(options)) {
            metadata(sce)[[key]] <- options[[key]]
        }
    }

    return(sce)
}

main <- function() {
    option_list <- list(
        make_option("--sample"),
        make_option("--dataset"),
        make_option("--platform", default="Visium"),
        make_option("--nPCs", type="integer", default=15),
        make_option("--logcounts", action="store_true", default=FALSE),
        make_option("--tenX", action="store_true", default=FALSE)
    )
    parser <- OptionParser(usage = "%prog [options] (counts rowData colData | --tenX dirname) RDS",
                           option_list=option_list)
    opt <- parse_args(parser, positional_arguments=c(2, 4))
    args <- opt$args
    options <- opt$options

    if (options$tenX) {
        sce <- make_SCE_from_10X(args[[1]])
        fout <- args[[2]]
    } else {
        sce <- make_SCE(args[[1]], args[[2]], args[[3]], options)
        fout <- args[[4]]
    }
    sce <- add_metadata(sce, options)
    sce <- spatialPreprocess(sce, platform=options$platform, n.PCs=options$nPCs)

    saveRDS(sce, fout)
}

if (sys.nframe() == 0) {
    main()
}
