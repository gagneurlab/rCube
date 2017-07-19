#' @title rCubeRate object and constructors
#' 
#' @description \code{rCubeRate} is a subclass of 
#' \code{RangedSummarizedExperiment}, used to store estimation results such as
#' synthesis and degradation rates
#' 
#' @rdname rCubeRate
setClass("rCubeRates", contains="RangedSummarizedExperiment")
