#' This function provides a wrapper function to estimate gene-specific
#' dispersion values.
#' 
#' @param experiment A \code{\link{rCubeExperiment}} object containing feature 
#' counts for labeled and total RNA-seq samples
#' @param method Either "DESeqDispMAP", "DESeqDispFit", "DESeqDispGeneEst", or 
#' "Replicate" for the type of fitting of dispersions to be used. For the three 
#' DESeq methods see als \code{\link{DESeq2}}
#' 
#' @seealso \code{\link{DESeq2}}
#' @return Returns an updated \code{\link{rCubeExperiment}} object with dispersion
#' estimated included into rowData information
#' @export
#' @author Leonhard Wachutka, Carina Demel
#' @examples 
#' data(geneCounts)
#' geneCounts <- estimateSizeDispersions(geneCounts, method='DESeqDispMAP')
#' rowRanges(geneCounts)
#'
#' geneCounts <- estimateSizeDispersions(geneCounts, method='Replicate')
#' rowRanges(geneCounts)
estimateSizeDispersions <- function(experiment, method = c('DESeqDispMAP','DESeqDispFit','DESeqDispGeneEst','Replicate'))
{
	if(method == 'Replicate')
	{
		return(estimateSizeDispersions_replicate(experiment))
	}else{
	    return(estimateSizeDispersionsDESeq(experiment, method))
	}
}

estimateSizeDispersions_replicate <- function(experiment)
{
	rowRanges(experiment)$dispersion = 40
	return(experiment)
}


#' Estimates gene-specific dispersion estimates for total and labeled samples,
#' using \code{\link[DESeq2]{DESeq2}}.
#'
#' @param experiment A \code{\link{rCubeExperiment}} object containing feature 
#' counts for labeled and total RNA-seq samples
#' @param method Type of disperion from \code{\link[DESeq2]{DESeq2}}, default 'DispMap'.
#'
#' @return Returns an updated \code{rCubeExperiment} object with
#' dispersion estimates for labeled and total samples included in the rowData.
#' @seealso \code{\link{DESeq2}}
#' @author Carina Demel
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#'
#' @examples
#' 
estimateSizeDispersionsDESeq <- function(experiment, method = c('DESeqDispMAP','DESeqDispFit','DESeqDispGeneEst'))
{
    counts <- assay(experiment)
    colData <- colData(experiment)
    labeledSample <- colData$LT
    
    expDesignTotal <- subset(colData, labeledSample == 'T')
    message(paste(
        "For the Total gene dispersion, the following samples are used:", 
        paste(rownames(expDesignTotal), collapse = ", ")))  
    deseqMatrixTotal <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(
        counts[ ,labeledSample == 'T'], colData = expDesignTotal,
        design = formula(~ condition)))
    ddsTotal <- suppressMessages(DESeq2::DESeq(deseqMatrixTotal))
    
    expDesignLabeled <- subset(colData, labeledSample == 'L')
    message(paste(
        "For the Labeled gene dispersion, the following samples are used:",
        paste(rownames(expDesignLabeled), collapse = ", ")))  
    deseqMatrixLabeled <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(
        counts[ ,labeledSample == 'L'], colData = expDesignLabeled,
        design = formula(~ condition)))
    ddsLabeled <- suppressMessages(DESeq2::DESeq(deseqMatrixLabeled))
    
    suppressWarnings(
        if (method == 'DESeqDispGeneEst') {
            dispersionLabeled <- mcols(ddsLabeled)$dispGeneEst
            dispersionTotal <- mcols(ddsTotal)$dispGeneEst
        } else if ( method == 'DESeqDispFit') {
            dispersionLabeled <- mcols(ddsLabeled)$dispFit
            dispersionTotal <- mcols(ddsTotal)$dispFit
        } else {
            dispersionLabeled <- mcols(ddsLabeled)$dispersion
            dispersionTotal <- mcols(ddsTotal)$dispersion
        }
    )
    # in cases where no dispersion could be fitted because of missing values,
    # set dispersion to median gene dispersion
    if (sum(is.na(dispersionLabeled)) > 0) {
        message(paste("Some feature-specific dispersions could not be",
                      "estimated (because of 0 counts) and were set to",
                      "the median dispersion value of all features."))
    }
    dispersionLabeled[is.na(dispersionLabeled)] <- 
        median(dispersionLabeled, na.rm=TRUE)
    dispersionTotal[is.na(dispersionTotal)] <- 
        median(dispersionTotal, na.rm=TRUE)
    
    rowRanges(experiment)$dispersion_L <- dispersionLabeled
    rowRanges(experiment)$dispersion_T <- dispersionTotal
    return(experiment)
}