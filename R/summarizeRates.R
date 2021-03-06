#' Summarize rates for features
#' 
#' @description Estimated rates are averaged for features which overlap with a 
#' top-level feature.
#'
#' @param featureRates An \code{rCubeRates} object, see constructor function 
#'  \code{\link{createRResultCubeRates}}.
#' @param topLevelFeatures \code{GRanges} annotation
#' @param by grouping method, one of 'mean', 'median'
#'
#' @return Returns a \code{rCubeExperiment} object which contains averaged rates
#' for the specified \code{topLevelFeatures}.
#' @export
#' @import GenomicRanges
#' @seealso \code{\link[GenomicRanges]{findOverlaps}}
#'
#' @examples
#' data(exonRates)
#' rows <- rowRanges(exonRates)
#' topLevelFeatures <- reduce(rows)
#' topLevelFeaturesRates <- summarizeRates(exonRates, topLevelFeatures, by='mean')
summarizeRates <- function(featureRates, topLevelFeatures, by=c('mean', 'median'))
{
    rows <- rowRanges(featureRates)
    rates <- assay(featureRates)
    ov <- findOverlaps(rows, topLevelFeatures)
    mergedRates <- .createRResultCubeRatesTopLevelFeatures(featureRates,
                                                            topLevelFeatures)
    mergeRates <- function(subject, rates, ov, by){
        hits <- queryHits(ov)[subjectHits(ov) == subject]
        if(length(hits) > 1){
            return(apply(rates[hits, ], 2, match.fun(by), na.rm=TRUE))
        }else{
            return(rates[hits, ])
        }
    }
    res <- t(sapply(unique(subjectHits(ov)), mergeRates, rates=rates, ov=ov, 
                    by=by))
    assay(mergedRates)[unique(subjectHits(ov)), ] <- res
    
    return(mergedRates)
}


## helper function to extract same rCubeRates setup as in featureRates but modify rowRanges
.createRResultCubeRatesTopLevelFeatures <- function(featureRates, topLevelFeatures)
{
    dm <- colData(featureRates)
    rates <- matrix(NA, nrow=length(topLevelFeatures), ncol=nrow(dm))
    
    se <- SummarizedExperiment(assays=list("rates"=rates),
                               rowRanges=topLevelFeatures,
                               colData=dm)
    colnames(se) <- dm$sample
    se <- new("rCubeRates", se)
    return(se)
}
