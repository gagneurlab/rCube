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
estimateSizeFactors <- function(featureCounts, spikeinCounts, 
                        method=c('spikeinGLM','spikeinMean','jointModel')){
    if(method=="spikeinGLM"){
        spikeinCounts <- calculateNormalizationBySpikeinGLM(spikeinCounts)
    }else if(method=="spikeinMean"){
        # calculateNormalizationBySpikeinMean
    }else{
        # calculateNormalizationByJointModel
    }
    
    # depends on how Leo's normalization methods work, which objects are updated
    return(spikeinCounts)
}



#' Fit GLM to spike-in read counts to extract sequencing depth and
#' cross-contamination rate for all samples
#'
#' @param spikeinCounts rCubeExperiment Object containing read counts for spike-ins
#'
#' @importFrom MASS glm.nb
#' 
#' @return rCubeExperiment Object for spike-ins with updated information
#' @export
#' @author Carina Demel
#'
#' @examples
#' data(spikeins)
#' calculateNormalizationBySpikeinGLM(spikeins)
calculateNormalizationBySpikeinGLM <- function(spikeinCounts){
 
    spikein.dataframe <- function(spikeinCounts){
        counts <- assay(spikeinCounts)
        rowData <- rowData(spikeinCounts)
        colData <- colData(spikeinCounts)
        
        counts.vec <- as.vector(counts)
        spikeins <- rowData$gene_id
        spikein.lengths <- rowData$length
        spikein.labeling <- rowData$labeled
        samples <- rownames(colData)
        conditions.labeling <- colData$labeled.sample
        
        mat <- data.frame(spike = rep(spikeins, length(samples)), 
                          length = rep(spikein.lengths, length(samples)),
                          spikein.labeled = rep(spikein.labeling, length(samples)),
                          sample = rep(paste(samples), each = length(spikeins)),
                          sample.labeling = rep(conditions.labeling, each = length(spikeins)),
                          counts = counts.vec)
        
        # additional columns: control for crosscontamination, with one value per
        # labeled sample and one value for ALL total samples (e.g.FALSE)
        # and log.length: natural logarithm of spike-in length
        mat$control.for.cross.contamination <- 
            ifelse(mat$sample.labeling == TRUE & mat$spikein.labeled == FALSE, TRUE, FALSE)
        mat$control.for.cross.contamination <- 
            factor(mat$control.for.cross.contamination)
        mat$sample.labeling <- factor(mat$sample.labeling)
        mat$ccc <- paste("L", rep(1:length(samples), each = length(spikeins)), 
                         collape = " ")
        mat$ccc[mat$control.for.cross.contamination == FALSE] <- "F"
        mat$ccc <- factor(mat$ccc)
        mat$log.length <- log(mat$length) #use natural logarithm
        
        return(mat)
    }
    mat <- spikein.dataframe(spikeinCounts)
    
    nb_model <- MASS::glm.nb(counts ~ offset(log.length) + sample + ccc + spike,
                             data=mat, link="log")
    fitted.counts <- nb_model$fitted.values

    coefs <- coefficients(nb_model)
    intercept <- coefs[1]
    
    seq.depths <- exp(coefs[grep("sample", names(coefs))]) 
    names(seq.depths) <- sub("sample", "", names(seq.depths))
    seq.depths <- c(1, seq.depths) #reference sample is 0, exp(0)=1
    names(seq.depths)[1] <- setdiff(samples, names(seq.depths))
    seq.depths <- seq.depths[samples]
    
    cross.cont.L <- exp(coefs[grep("cccL", names(coefs))])
    cross.cont <- rep(0, length(samples))
    names(cross.cont) <- samples
    cross.cont[conditions.labeling == "L"] <- cross.cont.L
    cross.cont[conditions.labeling == "T"] <- 1
    cross.cont <- cross.cont[samples]
    
    spikein.specific.bias <- c(1,exp(coefs[grep("spike",names(coefs))]))
    names(spikein.specific.bias) <- sort(spikeins)
    
    rowData(spikeinCounts)$spikein.specific.bias <- spikein.specific.bias
    rowData(spikeinCounts)$intercept <- intercept
    colData(spikeinCounts)[["sequencing.depth"]] <- seq.depths
    colData(spikeinCounts)[["cross-contamination"]] <- cross.cont
    assays(spikeinCounts)[["fitted.counts"]] <- matrix(fitted.counts, nrow=length(spikein.specific.bias))
    
    return(spikeinCounts)
}