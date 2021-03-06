#' @title rCubeExperiment object and constructors
#' 
#' @description \code{setupExperimentSpikeins} generates a \code{rCubeExperiment}
#' object for spike-ins, which can be used for counting and which is necessary
#' for normalization and estimation steps.
#'
#' @param rows GRanges containing spike-in/gene/exon/... annotation
#' @param designMatrix an optional data.frame containing sample information.
#' Either \code{designMatrix} or \code{files} needs to be present.
#' @param files An optional vector of sample file names in the format 
#' 'condition_L|T_labelingTime_Replicate.bam'. Either \code{designMatrix}
#' or \code{files} needs to be present.
#' @param length numeric vector with the length of each spike-ins
#' @param labelingState factor vector containing L and U indicating if a spikein
#' was labeled or unlabeled
#' @param counts an optional matrix with spike-in read counts
#'
#' @return An (empty) \code{rCubeExperiment} container
#' @export
#' @import SummarizedExperiment
#' @importFrom methods new
#' @importFrom stringr str_extract
#' @author Carina Demel, Leonhard Wachutka
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @rdname rCubeExperiment
#' 
#' @examples
#' data(spikeins)
#' data(spikeinLabeling)
#' spikeinLengths <- width(spikeins)
#' data(designMatrix)
#' spikeinCounts <- setupExperimentSpikeins(rows=spikeins, 
#' designMatrix=designMatrix, length=spikeinLengths, 
#' labelingState=spikeinLabeling)
setupExperimentSpikeins <- function(rows, designMatrix=NULL, files=NULL, 
                                    length=NULL, labelingState, counts=NULL){
    
    stopifnot(xor(is.null(designMatrix), is.null(files)) & !is.null(rows))
	
    if(is.null(designMatrix))
    {
        designMatrix <- .createDesignMatrix(files)
    }
    if(!is.null(length))
    {
        if(!is.null(names(length))){
            rows$length <- length[names(rows)]
        }else{
            rows$length <- length
        }
    }else{
        rows$length <- width(rows)
    }
    if(is.null(counts))
    {
        counts <- matrix(NA, nrow=length(rows), ncol=nrow(designMatrix))
        rownames(counts) <- names(rows)
    }
    if(!is.null(names(labelingState))){
        labelingState <- labelingState[names(rows)]
    }else{
        message("Your provided labelingState does not have names, rCube assumes that they are ordered according to the provided rows.")
    }
    
    
    rows$labelingState <- factor(labelingState)
    rows$labeledSpikein <- ifelse(rows$labelingState == "L", TRUE, FALSE)
    # rows$labeledSpikein <- factor(labelingState[names(rows)])
    
    
    rowData <- data.frame(length=rows$length, labelingState=labelingState, 
                            row.names=names(rows))
	
	if(is.null(names(rows)))
	{
		spikeins.SE <- SummarizedExperiment(assays=list("counts"=counts), rowRanges=rows, colData=designMatrix)
	}else{
		spikeins.SE <- SummarizedExperiment(assays=list("counts"=counts[names(rows), ]), rowRanges=rows, colData=designMatrix)
		
	}
    colnames(spikeins.SE) <- designMatrix$sample
    spikeins <- new("rCubeExperiment", spikeins.SE)
    return(spikeins)
}


#' @description \code{setupExperiment} generates an empty \code{rCubeExperiment}
#' container.
#' @export
#' @author Leonhard Wachutka
#' @rdname rCubeExperiment
#'
#' @examples
#' data(junctions)
#' bamfiles <- list.files(system.file("extdata/TimeSeriesExample/", package='rCube'), 
#' pattern=".bam$", full.names=TRUE)
#' junctionCounts <- setupExperiment(junctions, files=bamfiles)
setupExperiment <- function(rows, designMatrix=NULL, files=NULL){
    
    stopifnot(!(is.null(designMatrix) & is.null(files)))
    if(is.null(designMatrix))
    {
        designMatrix <- .createDesignMatrix(files)
    }
    
    counts <- matrix(NA, nrow=length(rows), ncol=nrow(designMatrix))
    se <- SummarizedExperiment(assays=list("counts"=counts), rowRanges=rows, colData=designMatrix)
    colnames(se) <- designMatrix$sample
    se <- new("rCubeExperiment", se)
    return(se)
}


.createDesignMatrix <- function(files)
{
    designMatrix <- data.table(filename=files)
    designMatrix[, stringr::str_extract(basename(filename), '([^_]+)')]
    designMatrix[, sample := stringr::str_extract(basename(filename), '([^_]+_){3}[^_^.]+')]
    designMatrix[, condition := as.factor(stringr::str_extract(sample, '[^_]+'))]
    designMatrix[, LT := as.factor(stringr::str_extract(sample, '(?<=_)[LT]+'))]
    designMatrix[, labelingTime := as.numeric(stringr::str_extract(sample, '(?<=_)[0-9]+'))]
    designMatrix[, replicate := as.factor(stringr::str_extract(sample, '(?<=_)[^_]*$'))]
    return(designMatrix)
}