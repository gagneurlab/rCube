#' Summarize rates for features
#' 
#' @description 
#'
#' @param featureRates An \code{rCubeRate object}, see constructor function 
#'  \code{\link{createRResultCubeRates}}.
#' @param topLevelFeatures \code{GRanges} annotation
#' @param method grouping method ## was this what we intended to do?
#'
#' @return Returns a \code{rCubeExperiment} object which contains summarized rates
#' for the specified topLevelFeatures.
#' @export
#'
#' @examples
#' data(geneRates)
#' rows <- rowRanges(geneRates)
#' topLevelFeatures <- reduce(rows)
#' topLevelFeaturesRates <- summarizeRates(geneRates, topLevelFeatures)
summarizeRates <- function(featureRates, topLevelFeatures, by=c('mean','median'))
{
    rows <- rowRanges(featureRates)
    rates <- assay(featureRates)
    ov <- findOverlaps(rows, topLevelFeatures)
    
    mergedRates <- .createRResultCubeRatesTopLevelFeatures(featureRates, topLevelFeatures, replicate)

    mergeRates <- function(subject, rates, ov, by){
        hits <- queryHits(ov)[subjectHits(ov) == subject]
        return(apply(rates, 2, mean, na.rm=T))#match.fun(by))) then na.rm does not work
    }
    res <- t(sapply(unique(subjectHits(ov)), mergeRates, rates=rates, ov=ov, by=by))
    assay(mergedRates) <- res
    
    return(mergedRates)
}


## helper function to extract same rCubeRates setup as in featureRates but modify rowRanges
.createRResultCubeRatesTopLevelFeatures <- function(featureRates, topLevelFeatures, replicate)
{
    dm = colData(featureRates)
    rates <- matrix(NA, nrow=length(topLevelFeatures), ncol=nrow(dm))
    
    se <- SummarizedExperiment(assays=list("rates"=rates), rowRanges=topLevelFeatures, colData=dm)
    colnames(se) = dm$sample
    se <- new("rCubeRates", se)
    return(se)
}