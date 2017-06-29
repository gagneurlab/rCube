#' Generate empty rCubeExperiment Object for spike-ins, which can be used for counting
#'
#' @param rows GRanges containing spike-in annotation
#' @param path string 
#' @param length numeric vector with the length of each spike-ins
#' @param labelingState factor vector or logical###??
#'
#' @return An empty rCubeExperiment container
#' @export
#' @import SummarizedExperiment
#' @author Carina Demel
#'
#' 
#' @examples
#' data(spikein.labeling)
#' data(spikein.lengths)
#' spikein.counts <- setupExperimentSpikeins(rows = spikeins, path = "", 
#' length = spikein.lengths, labelingState = spikein.labeling)
setupExperimentSpikeins <- function(rows, path, length, labelingState){
    rows$length <- length[names(rows)]
    rows$labelingState <- factor(labelingState[names(rows)]) ### TODO what to use here
    rows$labeled <- ifelse(rows$labelingState=="L", TRUE, FALSE) ### I use this now for Normalization.R
    rows$path <- path
    nfiles = length(list.files(path, pattern = ".bam"))
    #TODO ncol has to be set to final value directly
    counts <- matrix(NA, nrow = length(rows), ncol = nfiles)
    #todo what are the counts in the beginning?
    spikeins.SE <- SummarizedExperiment(assays = list("counts"=counts), rowRanges = rows)
    spikeins <- new("rCubeExperiment", spikeins.SE)
    return(spikeins)
}

#' Generate empty rCubeExperiment Object
#'
#' @param rows GRanges
#' @param designMatrix dgM 
#' 
#' @return An empty rCubeExperiment container
#' @export
#' @author Leonhard Wachutka
#'
#' @examples
#' 
setupExperiment <- function(rows, designMatrix){
	
	counts <- matrix(NA, nrow = length(rows), ncol = nrow(designMatrix))
	se <- SummarizedExperiment(assays = list("counts"=counts), rowRanges = rows,colData = designMatrix)
	colnames(se) = designMatrix$sample
	se <- new("rCubeExperiment", se)
	return(se)
}

