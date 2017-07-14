
#' 


#' @title Estimation of RNA rates
#' 
#' @description Estimation of synthesis and degradation rates (half-lives) based on
#' 4sU-labeled and total RNA-seq read counts.
#' 
#' @param featureCounts A \code{rCubeExperiment} object with feature count table
#' @param replicate list of character/factor vectors for which replicate or 
#' combination of replicates the results should be computed. If \code{NULL},
#' estimates for each replicate individually and all combinations of replicates
#' will be calculated.
#' @param method Type of estimation to be used. 'series' uses a set of TT-seq and
#' RNA-seq data sets with multiple labeling time points, 'single' works on individual
#' time points.
#' @param BPPARAM An instance of a \code{BiocParallelParam} class, e.g., 
#' \code{\link{MulticoreParam}}, \code{\link{SnowParam}}, \code{\link{DoparParam}}.
#' 
#' @return Returns a \code{rCubeRates} object with estimated synthesis and degradation
#' rates for each feature and sample and the specified replicate combinations
#' @export
#' @seealso \code{\link{BiocParallelParam}}
#' @import BiocParallel
#' @import data.table
#' @author Leonhard Wachutka, Carina Demel
#' @examples
#' ## estimate sequencing depths and cross-contamination values from spike-ins
#' data(spikeinCounts)
#' data(geneCounts)
#' geneCounts <- estimateSizeFactors(geneCounts, spikeinCounts, method='spikeinGLM')
#' ## estimate Dispersions for all genes
#' geneCounts <- estimateSizeDispersions(geneCounts, method='DESeqDispMAP')
#' ## estimate synthesis and degradation rates for individual replicates and combination
#' rates <- estimateRateByFirstOrderKineticsSingle(geneCounts, spikeinCounts, method='single', BPPARAM=NULL)
estimateRateByFristOrderKinetics = function(featureCounts, replicate, method=c('series','single'), BPPARAM=NULL)
{
	#replicate = c(1,2,'1:2')
	if(method == 'series')
	{
		rRates <- createRResultCubeRates(featureCounts, replicate)
		return(estimateRateByFristOrderKineticsSeries(featureCounts, rRates, BPPARAM=BPPARAM))
	}else{
	    return(.estimateRateByFristOrderKineticsSingle(featureCounts, replicate, BPPARAM=BPPARAM))
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