#' Title
#'
#' @param featureCounts rCubeExperiment Object containing read counts for features (genes, junctions, ...)
#' @param spikeinCounts rCubeExperiment Object containing read counts for spike-ins
#' @param method string specifying which normalization method should be used,
#' one of ('spikeinGLM','spikeinMean','jointModel')
#'
#' @return rCubeExperiment Object with updated metadata information
#' @export
#' @author Carina Demel
#'
#' @examples
#' 
estimateSizeFactors <- function(featureCounts, spikeinCounts, 
                        method=c('spikeinGLM','spikeinMean','jointModel')){
    if(method == "spikeinGLM"){
        spikeinCounts <- calculateNormalizationBySpikeinGLM(spikeinCounts)
    }else if(method == "spikeinMean"){
		featureCounts <- calculateNormalizationByMean(featureCounts, spikeinCounts)
    }else{
        # calculateNormalizationByJointModel
    }
    
    #' @Leo: do we return feature counts or spike-in counts?
    
    return(featureCounts)
}



#' Fit GLM to spike-in read counts to extract sequencing depth and
#' cross-contamination rate for all samples
#'
#' @param spikeinCounts A \code{rCubeExperiment} object containing read counts,
#' length and #' labeling status for spike-ins, as well as sample information
#'
#' @importFrom MASS glm.nb
#' 
#' @return An \code{\link{rCubeExperiment}} object for spike-ins with 
#' updated information in \code{RowData}
#'
#' @author Carina Demel
#' @examples
#' data(spikeinCounts)
#' spikeinCounts <- calculateNormalizationBySpikeinGLM(spikeinCounts)
calculateNormalizationBySpikeinGLM <- function(spikeinCounts){
    samples <- rownames(colData(spikeinCounts))
    spikeins <- rowRanges(spikeinCounts)$gene_id
    spikeinDataframe <- function(spikeinCounts){
        counts <- assay(spikeinCounts)
        rowData <- rowRanges(spikeinCounts)
        colData <- colData(spikeinCounts)
        
        countsVector <- as.vector(counts)
        spikeins <- rowData$gene_id
        spikeinLengths <- rowData$length
        spikeinLabeling <- rowData$labeledSpikein
        samples <- rownames(colData)
        conditionsLabeling <- colData$LT
        
        mat <- data.frame(spike = rep(spikeins, length(samples)), 
                          length = rep(spikeinLengths, length(samples)),
                          spikein.labeled = rep(spikeinLabeling, length(samples)),
                          sample = rep(paste(samples), each = length(spikeins)),
                          sample.labeling = rep(conditionsLabeling, each = length(spikeins)),
                          counts = countsVector)
        
        ## additional columns: control for crosscontamination, with one value per
        ## labeled sample and one value for ALL total samples (e.g.FALSE)
        ## and log.length: natural logarithm of spike-in length
        mat$control.for.crosscontamination <- 
            ifelse(mat$sample.labeling == "L" & mat$spikein.labeled == FALSE, TRUE, FALSE)
        mat$control.for.crosscontamination <- 
            factor(mat$control.for.crosscontamination)
        mat$sample.labeling <- factor(mat$sample.labeling)
        mat$ccc <- paste("L", rep(1:length(samples), each = length(spikeins)), 
                         collape = " ")
        mat$ccc[mat$control.for.crosscontamination == FALSE] <- "F"
        mat$ccc <- factor(mat$ccc)
        mat$log.length <- log(mat$length)
        
        return(mat)
    }
    mat <- spikeinDataframe(spikeinCounts)
    
    nb_model <- MASS::glm.nb(counts ~ offset(log.length) + sample + ccc + spike,
                             data=mat, link="log")
    fittedCounts <- nb_model$fitted.values

    coefs <- coefficients(nb_model)
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
    
    spikeinSpecificBias <- c(1,exp(coefs[grep("spike",names(coefs))]))
    names(spikeinSpecificBias) <- sort(spikeins)
    
    rowData(spikeinCounts)$spikeinSpecificBias <- spikeinSpecificBias
    rowData(spikeinCounts)$intercept <- intercept
    colData(spikeinCounts)$sequencing.depth <- sequencingDepths
    colData(spikeinCounts)$cross.contamination <- crossContamination
    assays(spikeinCounts)$fittedCounts <- matrix(fittedCounts,
                                                  nrow=length(spikeinSpecificBias))
    
    return(spikeinCounts)
}


calculateNormalizationByMean <- function(featureCounts, spikeinCounts){
	labled = subset(spikeinCounts, labelingState == 'L')
	lc = assay(labled)
	colData(featureCounts)$sizeFactor = colMeans(lc / rowSums(lc))
	return(featureCounts)
}