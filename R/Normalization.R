#' @title Estimation of sample-specific normalization factors
#' 
#' @description \code{estimateSizeFactors} provides a wrapper function to
#' estimate sample-specific size-factors by three different types. The 
#' 'spikeinGLM' method fits a GLM to spike-in read counts to extract sequencing 
#' depth values and cross-contamination rates for all samples.
#' 
#' TODO: explain other methods 
#' 
#' Additional informations for the fit can be found in the metadata of the 
#' resulting \code{rCubeExperiment}.
#'
#' @param featureCounts A \code{rCubeExperiment} object containing read
#' counts for features (genes, exons, introns, junctions, ...), see constructor 
#' function \code{\link{setupExperiment}}.
#' @param spikeinCounts A \code{rCubeExperiment} object containing read 
#' counts, length and labeling status for spike-ins, as well as sample 
#' information, see constructor function \code{\link{setupExperimentSpikeins}}.
#' @param method string specifying which normalization method should be used,
#' one of \code{('spikeinGLM', 'spikeinMean', 'jointModel')}.
#'
#' @return An updated \code{rCubeExperiment} object for features with sequencing
#'  depth values included in \code{colData} and updated metadata information.
#' 
#' @rdname estimateSizeFactors
#' @importFrom MASS glm.nb
#' @export
#' @author Carina Demel, Leonhard Wachutka
#'
#' @examples
#' # fitting a GLM to spike-in read counts
#' data(geneCounts)
#' data(spikeinCounts)
#' geneCounts <- estimateSizeFactors(geneCounts, spikeinCounts, method="spikeinGLM")
estimateSizeFactors <- function(featureCounts,
                                spikeinCounts,
                                method=c('spikeinGLM', 'spikeinMean', 'jointModel')){
    if(method == "spikeinGLM"){
        featureCounts <- .calculateNormalizationBySpikeinGLM(featureCounts, 
                                                             spikeinCounts)
    }else if(method == "spikeinMean"){
        featureCounts <- .calculateNormalizationByMean(featureCounts,
                                                       spikeinCounts)
    }else{
        ## TODO joint Model?
    }
    
    return(featureCounts)
}




.calculateNormalizationBySpikeinGLM <- function(featureCounts, spikeinCounts){
    samples <- rownames(colData(spikeinCounts))
    spikeins <- rowRanges(spikeinCounts)$gene_id
    conditionsLabeling <- colData(spikeinCounts)$LT
    
    mat <- .spikeinDataframe(spikeinCounts)
    
    nb_model <- MASS::glm.nb(counts ~ offset(log.length) + sample + ccc + spike,
                             data=mat, link="log")
    fittedCounts <- nb_model$fitted.values
    
    coefs <- stats::coefficients(nb_model)
    intercept <- coefs[1]
    
    sequencingDepths <- exp(coefs[grep("sample", names(coefs))]) 
    names(sequencingDepths) <- sub("sample", "", names(sequencingDepths))
    sequencingDepths <- c(1, sequencingDepths) ## reference sample is 0, exp(0)=1
    names(sequencingDepths)[1] <- setdiff(samples, names(sequencingDepths))
    sequencingDepths <- sequencingDepths[samples]
    
    crossContaminationLabeled <- exp(coefs[grep("cccL", names(coefs))])
    crossContamination <- rep(0, length(samples))
    names(crossContamination) <- samples
    crossContamination[conditionsLabeling == "L"] <- crossContaminationLabeled
    crossContamination[conditionsLabeling == "T"] <- 1
    crossContamination <- crossContamination[samples]
    
    spikeinSpecificBias <- c(1, exp(coefs[grep("spike", names(coefs))]))
    names(spikeinSpecificBias) <- sort(spikeins)
    
    resultsList <- list("sequencing.depth"=sequencingDepths,
                       "cross.contamination"=crossContamination,
                       "spikeinSpecificBias"=spikeinSpecificBias,
                       "intercept"=intercept,
                       "fittedCounts"=matrix(fittedCounts,
                                             nrow=length(spikeinSpecificBias)))
    colData(featureCounts)$sequencing.depth <- sequencingDepths
    colData(featureCounts)$cross.contamination <- crossContamination
    metadata(featureCounts) <- resultsList
    return(featureCounts)
}


.calculateNormalizationByMean <- function(featureCounts, spikeinCounts){
    labeled <- subset(spikeinCounts, labelingState == 'L')
    lc <- assay(labeled)
    colData(featureCounts)$sizeFactor <- colMeans(lc / rowSums(lc))
    return(featureCounts)
}


## helper function
.spikeinDataframe <- function(spikeinCounts){
    counts <- assay(spikeinCounts)
    rowData <- rowRanges(spikeinCounts)
    colData <- colData(spikeinCounts)
    
    countsVector <- as.vector(counts)
    spikeins <- rowData$gene_id
    spikeinLengths <- rowData$length
    spikeinLabeling <- rowData$labeledSpikein
    samples <- rownames(colData)
    conditionsLabeling <- colData$LT
    
    mat <- data.frame(spike=rep(spikeins, length(samples)),
                      length=rep(spikeinLengths, length(samples)),
                      spikein.labeled=rep(spikeinLabeling, length(samples)),
                      sample=rep(paste(samples), each=length(spikeins)),
                      sample.labeling=rep(conditionsLabeling, each=length(spikeins)),
                      counts=countsVector)
    
    ## additional columns: control for crosscontamination, with one value per
    ## labeled sample and one value for ALL total samples (e.g.FALSE)
    ## and log.length: natural logarithm of spike-in length
    mat$control.for.crosscontamination <- 
        ifelse(mat$sample.labeling == "L" & mat$spikein.labeled == FALSE, TRUE, FALSE)
    mat$control.for.crosscontamination <- 
        factor(mat$control.for.crosscontamination)
    mat$sample.labeling <- factor(mat$sample.labeling)
    mat$ccc <- paste("L", rep(1:length(samples), each=length(spikeins)), collape=" ")
    mat$ccc[mat$control.for.crosscontamination == FALSE] <- "F"
    mat$ccc <- factor(mat$ccc)
    mat$log.length <- log(mat$length)
    
    return(mat)
}