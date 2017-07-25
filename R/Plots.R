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
#' data(geneRates)
#' plotResultsReplicates(geneRates)
plotResultsReplicates <- function(featureRates){
    # nplots <- expand.grid(cond=unique(colData(featureRates)$condition), rate=c('synthesis', 'degradation'))
    
    
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
#' @param condition bla
#' @param replicate bla
#'
#'
#' @export
#' @importFrom graphics abline mtext par plot
#' @rdname diagnosticPlots
#'
#' @examples
#' data(spikeinCounts)
#' data(geneCounts)
#' geneCounts <- estimateSizeFactors(geneCounts, spikeinCounts, method='spikeinGLM')
#' data(geneRates)
#' plotFittedCounts(geneCounts, geneRates, "A", "1")
plotFittedCounts <- function(featureCounts, featureRates, condition=NULL, replicate=NULL){
    res.batches <- unique(data.frame(cond=colData(featureRates)$condition, rep=colData(featureRates)$replicate))
    res.batches <- res.batches[which(res.batches$cond %in% condition & res.batches$rep %in% replicate),]
    
    crossCont <- colData(featureCounts)$cross.contamination
    seqDepths <- colData(featureCounts)$sequencing.depth

    op <- par(no.readonly = TRUE)
    par(mfrow=c(nrow(res.batches), 2))
    
    for(r in 1:nrow(res.batches)){
        row <- res.batches[r, ]
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
        
        trueCountsL <- assay(featureCounts[, featureCounts$condition == row$cond &
                                            featureCounts$replicate %in% rep &
                                            featureCounts$LT == "L"])
        trueCountsT <- assay(featureCounts[, featureCounts$condition == row$cond &
                                            featureCounts$replicate %in% rep &
                                            featureCounts$LT == "T"])
        
        suppressWarnings(plot(expCountsL, trueCountsL, log="xy", 
                              xlab="Expected Counts", ylab="Observed Counts",
                              main=paste(row$cond, as.character(row$rep), 
                                         collapse=" ")))
        mtext("Labeled Counts", col="grey")
        abline(0, 1, col="grey")
        suppressWarnings(plot(expCountsT, trueCountsT, log="xy", 
                              xlab="Expected Counts", ylab="Observed Counts",
                              main=paste(row$cond, as.character(row$rep), 
                                         collapse=" ")))
        mtext("Total Counts", col="grey")
        abline(0, 1, col="grey")
    }
    par(op) ### Reset to previous plotting settings
}
