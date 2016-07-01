#' \code{M3D}: Uses the M3D method to find DMRs in bisufite sequencing data
#'
#' @docType package
#' @name M3D
NULL
#> NULL


#' @useDynLib M3D
#' @importFrom Rcpp sourceCpp
#' @import BiocGenerics
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @importFrom SummarizedExperiment colData rowRanges
#' @import BiSeq
#' @importFrom parallel mcmapply detectCores
NULL