#' @description \code{rCubeExperiment} is a subclass of 
#' \code{RangedSummarizedExperiment}, used to store read counts, annotations, and
#' sample metadata.
#' @rdname rCubeExperiment
#' @export
setClass("rCubeExperiment",contains = "RangedSummarizedExperiment",
         prototype=list(rowRanges=GRanges())
         )
