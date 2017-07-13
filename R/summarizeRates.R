#' Title
#'
#' @param featureRates An \code{rCubeRate object}, see constructor function 
#'  \code{\link{createRResultCubeRates}}.
#' @param topLevelFeatures \code{GRanges} annotation
#' @param by grouping
#' @param featureCounts A \code{rCubeExperiment} object, see constructor function 
#'  \code{\link{setupExperiment}}.
#'
#' @return Returns a \code{rCubeExperiment} object which contains summarized rates
#' for the specified topLevelFeatures.
#' @export
#'
#' @examples
summarizeRates <- function(featureRates, topLevelFeatures, by, featureCounts=NULL){
    
}
## merge different exons/junctions based on findOverlaps