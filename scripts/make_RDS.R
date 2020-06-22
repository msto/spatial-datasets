#!/usr/bin/env Rscript

library(optparse)
library(Matrix)
# suppressMessages(library(scater))
suppressMessages(library(SingleCellExperiment))

make_SCE_from_10X <- function(dirname) {
    spatial_dir <- file.path(dirname, "spatial")
    matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")

    if (!dir.exists(matrix_dir)) stop(paste0("Matrix directory does not exist: ", matrix_dir))
    if (!dir.exists(spatial_dir)) stop(paste0("Spatial directory does not exist: ", spatial_dir))

    colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header=F)
    # colnames(colData) <- c("spot", "in_tissue", "x", "y", "image_x", "image_y")
    colnames(colData) <- c("spot", "in_tissue", "col", "row", "imagecol", "imagerow")
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
make_SCE <- function(counts, rowData, colData) {
    # Load cleaned CSVs
    rowData <- read.csv(rowData, stringsAsFactors=FALSE)
    colData <- read.csv(colData, stringsAsFactors=FALSE)
    counts <- as.matrix(read.csv(counts, row.names=1, check.names=F, stringsAsFactors=FALSE))

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
        rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
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
    for (key in c("sample", "dataset")) {
        if (key %in% names(options)) {
            metadata(sce)[[key]] <- options[[key]]
        }
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
    sce <- scater::runPCA(sce, ncomponents=15)

    # De-noised PCA
    # TODO: solve error
    #   Error in density.default(x, adjust = adjust, from = min(x), to = max(x)) :
    #   need at least 2 points to select a bandwidth automatically
    # dec <- scran::modelGeneVarByPoisson(sce)
    # top <- scran::getTopHVGs(dec, prop=0.1)
    # sce <- scran::denoisePCA(sce, technical=dec, subset.row=top)

    # Other dimensionality reduction
    # sce <- scater::runTSNE(sce, dimred="PCA")
    # sce <- scater::runUMAP(sce, dimred="PCA")
}

main <- function() {
    option_list <- list(
        make_option("--sample"),
        make_option("--dataset"),
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
        sce <- make_SCE(args[[1]], args[[2]], args[[3]])
        fout <- args[[4]]
    }
    sce <- add_metadata(sce, options)
    sce <- process_SCE(sce, options)

    saveRDS(sce, fout)
}

if (sys.nframe() == 0) {
    main()
}
