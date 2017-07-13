#' Generate \code{\link{rCubeExperiment}} object for spike-ins, which can be
#' used for counting and which is necessary for normalization and estimation steps.
#'
#' @param rows GRanges containing spike-in annotation
#' @param designMatrix an optional data.frame containing sample information, see 
#' constructor functions \code{\link{createDesignMatrix}}. Either \code{designMatrix}
#' or \code{files} needs to be present.
#' @param files An optional vector of sample file names in the format 
#' 'condition_L|T_labelingTime_Replicate.bam'. Either \code{designMatrix}
#' or \code{files} needs to be present
#' @param length numeric vector with the length of each spike-ins
#' @param labelingState factor vector containing L and U indicating if a spikein
#' was labeled or unlabeled
#' @param counts an optional matrix with spike-in read counts
#'
#' @return An (empty) rCubeExperiment container
#' @export
#' @import SummarizedExperiment
#' @import stringr
#' @author Carina Demel, Leonhard Wachutka
#' @seealso \code{\link{SummarizedExperiment}}
#' 
#' @examples
#' data(spikeins)
#' data(spikeinLabeling)
#' data(spikeinLengths)
#' data(designMatrix)
#' spikeinCounts <- setupExperimentSpikeins(rows=spikeins, designMatrix=designMatrix, 
#' length=spikeinLengths, labelingState=spikein.labeling)
setupExperimentSpikeins <- function(rows, designMatrix = NULL, files = NULL, length = NULL, labelingState, counts=NULL){
	
	stopifnot(!(is.null(designMatrix) & is.null(files)))
	if(is.null(designMatrix))
	{
		designMatrix <- createDesignMatrix(files)
	}
	if(!is.null(length))
	{
		rows$length <- length[names(rows)]
	}else{
	    rows$length <- width(rows)
	}
    if(is.null(counts))
    {
        counts <- matrix(NA, nrow = length(rows), ncol = nrow(designMatrix))
        rownames(counts) <- names(rows)
    }
    
    rows$labelingState <- factor(labelingState[names(rows)])
    rows$labeledSpikein <- ifelse(rows$labelingState=="L", TRUE, FALSE) ### I use this now for Normalization.R
 
    rowData <- data.frame(length=rows$length, labelingState=labelingState, row.names=names(rows))
    spikeins.SE <- SummarizedExperiment(assays=list("counts"=counts[names(rows),]), rowRanges=rows, colData=designMatrix)
    colnames(spikeins.SE) <- designMatrix$sample
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
	
	stopifnot(!(is.null(designMatrix) & is.null(files)))
	if(is.null(designMatrix))
	{
		designMatrix <- createDesignMatrix(files)
	}
	
	counts <- matrix(NA, nrow=length(rows), ncol=nrow(designMatrix))
	se <- SummarizedExperiment(assays=list("counts"=counts), rowRanges=rows, colData=designMatrix)
	colnames(se) <- designMatrix$sample
	se <- new("rCubeExperiment", se)
	return(se)
}

createDesignMatrix = function(files)
{
	designMatrix <- data.table(filename=files)
	designMatrix[,str_extract(basename(filename),'([^_]+)')]
	designMatrix[,sample:=str_extract(basename(filename),'([^_]+_){3}[^_^.]+')]
	designMatrix[,condition:=as.factor(str_extract(sample,'[^_]+'))]
	designMatrix[,LT:=as.factor(str_extract(sample,'(?<=_)[LT]+'))]
	designMatrix[,labelingTime:=as.numeric(str_extract(sample,'(?<=_)[0-9]+'))]
	designMatrix[,replicate:=as.factor(str_extract(sample,'(?<=_)[^_]*$'))]
	return(designMatrix)
}