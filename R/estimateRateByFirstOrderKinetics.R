#' @title Estimation of RNA rates
#' 
#' @description Estimation of synthesis and degradation rates (half-lives) based
#' on 4sU-labeled and total RNA-seq read counts.
#' 
#' @param featureCounts A \code{rCubeExperiment} object with feature count table
#' @param replicate A vector of integers/strings indicating for which individual
#' replicates and for which combination of replicates the results should be 
#' calculated.
#' For specifying combinations of replicates, please separate them by a ':'.
#' If \code{NULL}, estimates for each replicate individually and all combinations 
#' of replicates will be calculated.
#' @param method Type of estimation to be used. 'single' works on individual
#' time points, 'series' uses a set of TT-seq and
#' RNA-seq data sets with multiple labeling time points.
#' @param BPPARAM An instance of a \code{BiocParallelParam} class, e.g., 
#' \code{\link{MulticoreParam}}, \code{\link{SnowParam}}, \code{\link{DoparParam}}.
#' 
#' @return Returns a \code{rCubeRates} object with estimated synthesis and 
#' degradation rates for each feature and sample and the specified replicate 
#' combinations.
#' 
#' @export
#' @seealso \code{\link{BiocParallelParam}}
#' @import BiocParallel
#' @importFrom data.table as.data.table rbindlist
#' @importFrom stats dnbinom formula optim runif
#' @importFrom utils as.relistable capture.output
#' 
#' @author Leonhard Wachutka, Carina Demel
#' @examples
#' ## estimate size factors and cross-contamination values from spike-ins
#' data(spikeinCounts)
#' data(exonCounts)
#' exonCounts <- estimateSizeFactors(exonCounts, spikeinCounts, method='spikeinGLM')
#' 
#' ## estimate dispersions for all genes
#' exonCounts <- estimateSizeDispersions(exonCounts, method='DESeqDispMAP')
#' 
#' ## set number of iterations for random initialization
#' elementMetadata(exonCounts)$numberOfInterations <- 1
#' 
#' ## estimate synthesis and degradation rates for individual replicates and combination
#' rates <- estimateRateByFirstOrderKinetics(exonCounts, replicate=c(1, 2, "1:2"),
#' method='single', BPPARAM=BiocParallel::MulticoreParam(1))
estimateRateByFirstOrderKinetics <- function(featureCounts,
                                             replicate,
                                             method=c('single', 'series'),
                                             BPPARAM=NULL)
{
    if(length(method) > 1){
        method <- method[1]
        message("More than one method was given and only the first element will be used.")
    }
    
    if(method == 'series')
    {
        rRates <- createRResultCubeRates(featureCounts, replicate)
        return( .estimateRateByFirstOrderKineticsSeries(featureCounts,
                                                        rRates,
                                                        BPPARAM=BPPARAM))
    }else{
        return( .estimateRateByFirstOrderKineticsSingle(featureCounts,
                                                        replicate,
                                                        BPPARAM=BPPARAM))
    }
}


#' @description \code{createRResultCubeRates} is a constructor for 
#' \code{rCubeExperiment} objects which contains slots for synthesis and 
#' degradation rates for all replicates/combination of replicates.
#'
#' @param featureCounts a \code{rCubeExperiment} object with set 'condition'
#' column in \code{colData}.
#' @param replicate a vector of integers and strings indicating single replicates/
#' combination of replicates for which the results should be calculated.
#'
#' @return Returns an empty rCubeRate object.
#' @export
#' @importFrom methods new
#' @author Leonhard Wachutka
#' @rdname rCubeRate
#'
#' @examples
#' data(exonCounts)
#' rates <- createRResultCubeRates(exonCounts, replicate=c(1, 2, '1:2'))
#' rates
createRResultCubeRates <- function(featureCounts, replicate)
{
    #designmatrix
    cond <- unique(featureCounts$condition)
    dm <- expand.grid(condition=cond, 
                      rate=as.factor(c('synthesis', 'degradation')), 
                      replicate=as.character(replicate))
    dm$sample <- paste0(dm$condition, '_', dm$rate, '_', dm$replicate)
    
    rates <- matrix(NA, nrow=length(featureCounts), ncol=nrow(dm))
    
    se <- SummarizedExperiment(assays=list("rates"=rates), 
                               rowRanges=rowRanges(featureCounts), colData=dm)
    colnames(se) <- dm$sample
    se <- new("rCubeRates", se)
    return(se)
}


#' @description \code{createRResultCubeRatesExtended} is a constructor for 
#' \code{rCubeExperiment} objects which contains slots for synthesis and 
#' degradation rates, as well as half.lives, labeled and unlabeled RNA amount 
#' estimates for all replicates/combination of replicates.
#'
#' @export
#' @rdname rCubeRate
#'
#' @examples
#' #' data(exonCounts)
#' rates <- createRResultCubeRatesExtended(exonCounts, replicate=c(1, 2, '1:2'))
#' rates
createRResultCubeRatesExtended <- function(featureCounts, replicate)
{
    #designmatrix
    cond <- unique(featureCounts$condition)
    dm <- expand.grid(condition=cond, rate=as.factor(c('synthesis', 'degradation', 'half.life', 'labeled.amount', 'unlabeled.amount')), replicate=as.character(replicate))
    dm$sample <- paste0(dm$condition, '_', dm$rate, '_', dm$replicate)
    
    rates <- matrix(NA, nrow=length(featureCounts), ncol=nrow(dm))
    
    se <- SummarizedExperiment(assays=list("rates"=rates), 
                               rowRanges=rowRanges(featureCounts), colData=dm)
    colnames(se) <- dm$sample
    se <- new("rCubeRates", se)
    return(se)
}
