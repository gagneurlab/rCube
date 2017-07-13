#' Title
#' 
#' @return A \code{rCubeRates} object with rates
#' @export
#' @author Leonhard Wachutka
#'
estimateRateByFristOrderKinetics = function(featureCounts, replicate, method=c('series','single'), BPPARAM=NULL)
{
	#replicate = c(1,2,'1:2')
	if(method == 'series')
	{
		rRates <- createRResultCubeRates(featureCounts, replicate)
		return(estimateRateByFristOrderKineticsSeries(featureCounts, rRates, BPPARAM=BPPARAM))
	}else{
	    return(estimateRateByFristOrderKineticsSingle(featureCounts, replicate, BPPARAM=BPPARAM))
	}
}


#' @description \code{createRResultCubeRates} is a constructor for \code{rCubeExperiment}
#' objects which contains slots for synthesis and degradation rates for all
#' replicates/combination of replicates.
#'
#' @param featureCounts a \code{rCubeExperiment} object with set 'condition'
#' column in \code{colData}.
#' @param replicate a list of single replicates/combination of replicates for
#' which the results should be calculated.
#'
#' @return Returns an empty rCubeRate object.
#' @export
#' @author Leonhard Wachutka
#' @rdname rCubeRate
#'
#' @examples
#' data(geneCounts)
#' rates <- createRResultCubeRates(geneCounts, list('1'=1,'2'=2,'1,2'=c(1,2)))
#' rates
createRResultCubeRates = function(featureCounts, replicate)
{
	#designmatrix
	cond <- unique(featureCounts$condition)
	dm <- expand.grid(condition=cond, rate=as.factor(c('synthesis', 'degradation')), replicate=as.character(replicate))
	dm$sample <- paste0(dm$condition, '_', dm$rate, '_', dm$replicate)
	
	rates <- matrix(NA, nrow=length(featureCounts), ncol=nrow(dm))
	
	se <- SummarizedExperiment(assays=list("rates"=rates), rowRanges=rowRanges(featureCounts), colData=dm)
	colnames(se) <- dm$sample
	se <- new("rCubeRates", se)
	return(se)
	
}

#' @description \code{createRResultCubeRatesExtended} is a constructor for 
#' \code{rCubeExperiment} objects which contains slots for synthesis and degradation rates,
#' as well as half.lives, labeled and unlabeled RNA amount estimates for all
#' replicates/combination of replicates. Used for \code{\link{estimateRateByFirstOrderKineticsSingle}}.
#'
#' @export
#' @rdname rCubeRate
#'
#' @examples
#' #' data(geneCounts)
#' rates <- createRResultCubeRatesExtended(geneCounts, list('1'=1,'2'=2,'1,2'=c(1,2)))
#' rates
createRResultCubeRatesExtended <- function(featureCounts, replicate)
{
    #designmatrix
    cond <- unique(featureCounts$condition)
    dm <- expand.grid(condition=cond, rate=as.factor(c('synthesis', 'degradation', 'half.life', 'labeled.amount', 'unlabeled.amount')), replicate=as.character(replicate))
    dm$sample <- paste0(dm$condition, '_', dm$rate, '_', dm$replicate)
    
    rates <- matrix(NA, nrow=length(featureCounts), ncol=nrow(dm))
    
    se <- SummarizedExperiment(assays=list("rates"=rates), rowRanges=rowRanges(featureCounts), colData=dm)
    colnames(se) <- dm$sample
    se <- new("rCubeRates", se)
    return(se)
    
}