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
#' @import stringr
#' @author Carina Demel, Leonhard Wachutka
#'
#' 
#' @examples
#' data(spikein.labeling)
#' data(spikein.lengths)
#' spikein.counts <- setupExperimentSpikeins(rows = spikeins, path = "", 
#' length = spikein.lengths, labelingState = spikein.labeling)
setupExperimentSpikeins <- function(rows, designMatrix = NULL, files = NULL, length = NULL, labelingState){
	
	stopifnot(!(is.null(designMatrix)&is.null(files)))
	if(is.null(designMatrix))
	{
		designMatrix = createDesignMatrix(files)
	}
	if(!is.null(length))
	{
		rows$length <- length[names(rows)]
	}
    
    rows$labelingState <- factor(labelingState) ### TODO what to use here
    rows$labeled <- ifelse(rows$labelingState=="L", TRUE, FALSE) ### I use this now for Normalization.R
    #rows$path <- path
    
    #TODO ncol has to be set to final value directly
    counts <- matrix(NA, nrow = length(rows), ncol = nrow(designMatrix))
    #todo what are the counts in the beginning?
    spikeins.SE <- SummarizedExperiment(assays = list("counts"=counts), rowRanges = rows, colData = designMatrix)

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
setupExperiment <- function(rows, designMatrix = NULL, files = NULL){
	
	stopifnot(!(is.null(designMatrix)&is.null(files)))
	if(is.null(designMatrix))
	{
		designMatrix = createDesignMatrix(files)
	}
	
	counts <- matrix(NA, nrow = length(rows), ncol = nrow(designMatrix))
	se <- SummarizedExperiment(assays = list("counts"=counts), rowRanges = rows,colData = designMatrix)
	colnames(se) = designMatrix$sample
	se <- new("rCubeExperiment", se)
	return(se)
}

createDesignMatrix = function(files)
{
	designMatrix = data.table(filename = files)
	designMatrix[,str_extract(basename(filename),'([^_]+)')]
	designMatrix[,sample:=str_extract(basename(filename),'([^_]+_){3}[^_^.]+')]
	designMatrix[,condition := as.factor(str_extract(sample,'[^_]+'))]
	designMatrix[,LT := as.factor(str_extract(sample,'(?<=_)[LT]+'))]
	designMatrix[,labelingTime := as.numeric(str_extract(sample,'(?<=_)[0-9]+'))]
	designMatrix[,replicate:=as.factor(str_extract(sample,'(?<=_)[^_]*$'))]
	return(designMatrix)
}