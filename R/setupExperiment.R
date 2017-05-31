#' Generate empty rCubeExperiment Object for spike-ins, which can be used for counting
#'
#' @param rows GRanges containing spike-in annotation
#' @param path string 
#' @param length numeric vector with the length of each spike-ins
#' @param labelingState factor vector
#'
#' @return An empty rCubeExperiment container
#' @export
#' @author Carina Demel
#'
#' @examples
#' data(spikein.labeling)
#' data(spikein.lengths)
#' spikein.counts <- setupExperimentSpikeins(rows = spikeins, path = "", 
#' length = spikein.lengths, labelingState = spikein.labeling)
setupExperimentSpikeins <- function(rows, path, length, labelingState){
    rows$length <- length[names(rows)]
    rows$labelingState <- factor(labelingState[names(rows)])
    rows$path <- path
    counts <- matrix(NA, nrow = length(rows), ncol = 1)
    #todo what are the counts in the beginning?
    spikeins.SE <- SummarizedExperiment(assays = counts, rowRanges = rows)
    spikeins <- new("rCubeExperiment", spikeins.SE)
    return(spikeins)
}