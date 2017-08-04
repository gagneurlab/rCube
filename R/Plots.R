# TODO Plot: model time vs normalization factor

#' @title Diagnostic plots 
#'
#' @description \code{plotSpikeinCountsVsSample} visualizes spike-in read counts
#' per sample to check if experiment worked and labeled spike-ins are enriched,
#' but unlabeled spike-ins are depleted in labeled samples, whereas all spike-ins
#' should be present to a similar extend in total RNA-seq samples.
#'
#' @param spikeinCounts A \code{rCubeExperiment} object containing read 
#' counts and labeling status for spike-ins, as well as sample information,
#' see constructor function \code{\link{setupExperimentSpikeins}}.
#'
#' @return A \code{trellis} object.
#' @importFrom data.table melt
#' @importFrom grDevices colorRampPalette
#' @import lattice
#' @export
#' @author Carina Demel
#' @rdname diagnosticPlots
#'
#' @examples
#' data(spikeinCounts)
#' plotSpikeinCountsVsSample(spikeinCounts)
plotSpikeinCountsVsSample <- function(spikeinCounts)
{
    counts <- assay(spikeinCounts)
    sample <- colData(spikeinCounts)$sample
    rows <- rowRanges(spikeinCounts)
    spikein.names <- names(rows)
    rownames(counts) <- spikein.names
    
    df <- melt(counts, value.name='counts')
    
    labeledPalette <- colorRampPalette(c("goldenrod1", "orange", "red"), space="rgb")
    unlabeledPalette <- colorRampPalette(c("darkslategray1", "blue"), space="rgb")
    spikein.colors <- c(labeledPalette(sum(rows$labeledSpikein == TRUE)), 
                       unlabeledPalette(sum(rows$labeledSpikein == FALSE)))
    
    names(spikein.colors) <- c(spikein.names[rows$labeledSpikein == TRUE], 
                               spikein.names[rows$labeledSpikein == FALSE])
    xyplot(counts ~ Var2, scales=list(y=list(log=2)), data=df, pch=16, 
           col=spikein.colors[as.character(df$Var1)], ylab="Spike-in Counts", 
           xlab="Sample", key=list(space="bottom", columns=2,
                                   points=list(col=spikein.colors, pch=16),
                                   text=list(names(spikein.colors))
           ))
}



#' @description \code{plotResultsReplicates} 
#'
#'
#' @export
#' @rdname diagnosticPlots
#' @import graphics
#'
#' @examples
#' data(exonRates)
#' plotResultsReplicates(exonRates)
plotResultsReplicates <- function(featureRates){
    # nplots <- expand.grid(cond=unique(colData(featureRates)$condition), rate=c('synthesis', 'degradation'))
    
    #TODO
    
    for(cond in unique(colData(featureRates)$condition)){
        for(rate in c('synthesis', 'degradation')){
            rates <- assay(featureRates)[, featureRates$rate == rate & featureRates$condition == cond]
            pairs(rates)
        }
    }
    
}


#' @description \code{plotFittedCounts} plots for all investigated conditions and
#' replicates the observed vs expected read counts in labeled and total RNA-seq
#' samples.
#'
#' @param featureCounts A \code{rCubeExperiment} object containing read counts
#' and sample information.
#' @param featureRates A \code{rCubeRates} object containing estimated labeled
#' and total RNA amounts for all desired replicates/conditions.
#' @param condition (optional) Character vector indicating subset of samples
#' for which plot should be produced
#' @param replicate (optional) replicate/combination for which plot should be produced
#'
#'
#' @export
#' @import ggplot2
#' @rdname diagnosticPlots
#'
#' @examples
#' data(exonRates)
#' plotFittedCounts(exonCounts, exonRates, "ActivatedJurkat", "1")
plotFittedCounts <- function(featureCounts, featureRates, condition=NULL, replicate=NULL){
    res.batches <- unique(data.frame(cond=colData(featureRates)$condition, rep=colData(featureRates)$replicate))
    
    if(is.null(condition)){
        condition <- res.batches$cond
    }
    if(is.null(replicate)){
        replicate <- res.batches$rep
    }
    res.batches <- res.batches[which(res.batches$cond %in% condition & res.batches$rep %in% replicate),]
    
    crossCont <- colData(featureCounts)$cross.contamination
    seqDepths <- colData(featureCounts)$sequencing.depth
    
    reps = do.call("rbind", sapply(1:nrow(res.batches), function(r){ rep=unlist(strsplit(as.character(res.batches[r,]$rep), ':'));
    return(suppressWarnings(data.frame(res.batches[r,], repind=rep)))
    }, simplify = FALSE))
    
    # counts <- expand.grid(sample=paste(res.batches$cond, res.batches$rep), feature = 1:nrow(rowData(featureCounts)))
    counts <- expand.grid(sample=paste(reps$cond, reps$rep), feature=1:nrow(rowData(featureCounts)), LT=c("L","T"))
    # counts <- cbind(counts, trueCountsL=NA, trueCountsT=NA, expCountsL=NA, expCountsT=NA)
    counts <- cbind(counts, trueCounts=NA, expCounts=NA)
    
    for(r in 1:nrow(res.batches)){
        row <- res.batches[r, ]
        p <- paste(row$cond, row$rep)
        labeledAmount <- assay(featureRates[, featureRates$condition == row$cond & featureRates$replicate == row$rep & featureRates$rate == 'labeled.amount'])
        unlabeledAmount <- assay(featureRates[, featureRates$condition == row$cond & featureRates$replicate == row$rep & featureRates$rate == 'unlabeled.amount'])
        
        rep <- unlist(strsplit(as.character(row$rep), ':'))
        labeledSamples <- which(colData(featureCounts)$condition == row$cond & colData(featureCounts)$replicate %in% rep & colData(featureCounts)$LT == "L")
        totalSamples <- which(colData(featureCounts)$condition == row$cond & colData(featureCounts)$replicate %in% rep & colData(featureCounts)$LT == "T")
        
        expCountsT <- .vgetExpectedCounts(width(featureCounts), labeledAmount,
                                          unlabeledAmount, N=1,
                                          crossCont[totalSamples],
                                          seqDepths[totalSamples])
        expCountsL <- .vgetExpectedCounts(width(featureCounts), labeledAmount,
                                          unlabeledAmount, N=1,
                                          crossCont[labeledSamples],
                                          seqDepths[labeledSamples])
        
        counts[counts$sample == p & counts$LT == "L", "expCounts"] <- as.vector(t(expCountsL))
        counts[counts$sample == p & counts$LT == "T", "expCounts"] <- as.vector(t(expCountsT))
        
        trueCountsL <- assay(featureCounts[, featureCounts$condition == row$cond &
                                               featureCounts$replicate %in% rep &
                                               featureCounts$LT == "L"])
        trueCountsT <- assay(featureCounts[, featureCounts$condition == row$cond &
                                               featureCounts$replicate %in% rep &
                                               featureCounts$LT == "T"])
        counts[counts$sample == p & counts$LT == "L", "trueCounts"] <- as.vector(t(trueCountsL))
        counts[counts$sample == p & counts$LT == "T", "trueCounts"] <- as.vector(t(trueCountsT))
    }
    suppressWarnings(print(ggplot(counts, aes(x=expCounts, y=trueCounts)) + 
                               geom_point() + 
                               facet_grid(LT ~ sample) + 
                               scale_x_log10() + 
                               scale_y_log10() + 
                               geom_abline(intercept=0, slope=1) + 
                               xlab("Expected read counts") +
                               ylab("Observed read counts") +
                               ggtitle("Expected vs observed read counts")))
}
