#' @title Adjust synthesis rate to transcripts per cell per minute
#' 
#' @description Based on the knowledge how many spike-ins were originally added 
#' to a defined
#' number of cells, estimated synthesis rates can be adjusted to retrieve
#' a rate for transcripts per cell and minute.
#' 
#' @param spikeinCounts \code{rCubeExperiment} object with spikeins read counts
#' @param spikeinGenome A \code{DNAStringSet} containing the spike-in sequences.
#' Names must be identical to the \code{seqnames} in the \code{spikeinCounts} object.
#' @param featureRates \code{rCubeRates} object with resulting synthesis rates 
#' that should be adjusted.
#' @param ratio_4sU Percentage of 4sU vs U volume for reverse transcription of 
#' spike-ins.
#' @param spikeinAmount amount of spike-ins added to each sample in ng. Either 
#' a single value, if same amounts were used for all spike-ins, or a vector.
#' Alternatively, a column in the \code{rowData} of the \code{spikeinCounts} can
#' be set, named 'spikeinAmount'.
#' @param cellNumber Number of cells used for the experiment. Either a single 
#' value if same cell numbers were used for all samples, or a vector.
#' Alternatively, a column in the \code{colData} of the \code{spikeinCounts} can
#' be set, named 'cellNumber'.
#' Note: if you want to use different spike-ins amounts for different samples
#' with different cell numbers,
#' please split the data according to the samples and run it separately.
#' 
#' @return Returns a matrix of synthesis rates, normalized to 
#' transcripts/cell/minute.
#' 
#' @author Björn Schwalb, Carina Demel
#' @export
#' @import Biostrings
#' @examples 
#' data("spikeinGenome")
#' data("spikeinCounts")
#' # due to downsampling of read counts to 5% of originial counts, we need to adjust here:
#' rowRanges(spikeinCounts)$spikeinAmount <- 25*0.05
#' colData(spikeinCounts)$cellNumber <- 5*10^7
#' data("exonRates")
#' adjustedSynthesis <- adjustSynthesisRate(spikeinCounts, spikeinGenome, exonRates)
adjustSynthesisRate <- function(spikeinCounts, spikeinGenome, featureRates, ratio_4sU=0.1, spikeinAmount=1, cellNumber=10^7){
    #molecular weights [g/mol]:
    Mw_A <- 329.2
    Mw_U <- 306.2
    Mw_4sU <- 322.26 #U + S - O: 306.2 + 32.06 - 16 = 322.26
    Mw_C <- 305.2
    Mw_G <- 345.2
    Mw_triposphate <- 159 # 5' triphosphate
    
    if(!is.null(metadata(spikeinCounts)$ratio_4sU)){
        ratio_4sU <- metadata(spikeinCounts)$ratio_4sU
    }
    stopifnot( (ratio_4sU <= 1) )
    
    rows <- rowRanges(spikeinCounts)
    if(!is.null(rows$spikeinAmount)){
        spikeinAmount <- rows$spikeinAmount
    }
    cols <- colData(spikeinCounts)
    if(!is.null(cols$cellNumber)){
        cellNumber <- cols$cellNumber
    }
    
    #reorder
    spikeinGenome <- spikeinGenome[match(seqnames(rows), names(spikeinGenome))]
    
    # count nucleotides in DNA, but T is U in RNA
    nucleotideCounts <- Biostrings::letterFrequency(spikeinGenome, c("A", "C", "G", "T"))
    rownames(nucleotideCounts) <- names(spikeinGenome)
    nucleotideCounts = cbind(nucleotideCounts, nucleotideCounts[,"T"]) #for 4sU
    spikesChrs <- rownames(nucleotideCounts)
    
    # for all spikeins
    ratio_4sU_spikeins = ifelse(rows$labeledSpikein == TRUE, ratio_4sU, 0)
    spikeinMolecularWeight <- rowSums(nucleotideCounts * cbind(Mw_A, Mw_C, Mw_G, (1-ratio_4sU_spikeins)*Mw_U, ratio_4sU_spikeins*Mw_4sU)) + Mw_triposphate
    
    # only labeled spikeins
    labeledSpikes = names(rows)[which(rows$labeledSpikein == TRUE)]
    labeledSpikesChrs = spikesChrs[which(rows$labeledSpikein == TRUE)]
    
    # Avogadro constant Na [1/mol] = N /n [mol]
    # => N = Na * n
    # Molar Mass M  = mass m [g] / n [mol]
    # => n = m[g] / M [g/mol]
    # ==> N = Na [1/mol] * m [g] / M  [g/mol]
    AvogadroConst <- 6.02214085774*10^23
    #nanogram to gram conversion:  * 10^(-9)
    spikeinNumber <- spikeinAmount * 10^(-9) * AvogadroConst / spikeinMolecularWeight
    spikeinNumber <- spikeinNumber[labeledSpikesChrs]
    
    # use counts only from labeled samples as synthesis is derived from labeled samples
    spikeinLength <- width(rows[labeledSpikes])
    counts <- assay(spikeinCounts)
    labeledSamples = colData(spikeinCounts)$LT == "L"
    if(length(cellNumber) == 1){
        cellNumber <- rep(cellNumber, length(labeledSamples))
    }
    labelingTime <- cols$labelingTime[labeledSamples]

    spikeinNumberPerCell <- outer(spikeinNumber, cellNumber[labeledSamples], FUN="/")
    # median over spike-ins
    conversionFactorToAmountPerCell = apply(t(t(counts[labeledSpikes, labeledSamples])/labelingTime) / (spikeinLength * spikeinNumberPerCell), 2, median)
    # mean over replicates for robustness, split by sample!
    conversionFactorToAmountPerCell = sapply(split(conversionFactorToAmountPerCell, cols$condition[labeledSamples]), mean)
    # message("Conversion Factor to amount per cell:")
    # message(paste(paste(names(conversionFactorToAmountPerCell),conversionFactorToAmountPerCell, sep=": ") , collapse=", "))

    #for each condition
    adjustedSR <- do.call("cbind", lapply(names(conversionFactorToAmountPerCell), function(cond){
        assay(featureRates[, which(featureRates$rate == 'synthesis' & featureRates$condition == cond)]) / conversionFactorToAmountPerCell[cond]
    }))
    return(adjustedSR)
}