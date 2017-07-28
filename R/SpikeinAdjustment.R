#' Title
#' TODO
#'
#' @param spikeinCounts \code{rCubeExperiment} object with spikeins read counts
#' @param featureRates \code{rCubeRates} object with resulting synthesis rates that
#' should be adjusted
#' @param ratio_4sU Percentage of 4sU vs U volume for reverse transcription of spike-ins
#'
#' @author Carina Demel
#' @importFrom utils data
adjustSynthesisRate <- function(spikeinCounts, featureRates, ratio_4sU=0.1){
    #molecular weights:
    Mw_A = 329.2
    Mw_U = 306.2
    Mw_4sU = 322.26 #U + S - O: 306.2 + 32.06 - 16 = 322.26
    Mw_C = 305.2
    Mw_G = 345.2
    Mw_triposphate = 159 # 5' triphosphate
    
    stopifnot( (ratio_4sU <= 1) )
    
    if(!is.null(metadata(spikeinCounts)$ratio_4sU)){
        ratio_4sU = metadata(spikeinCounts)$ratio_4sU
    }
    data("spikeinGenome", envir=environment())
    nucleotide.counts = letterFrequency(spikein.genome, c("A", "C", "G", "T"))
    # spikein.molecular.weight = c(nucleotide.counts[c("chrS2","chrS4","chrS8"),"A"] * Mw_A +
    #                                  (1-ratio_4sU) * nucleotide.counts[c("chrS2","chrS4","chrS8"),"T"] * Mw_U +
    #                                  nucleotide.counts[c("chrS2","chrS4","chrS8"),"C"] * Mw_C +
    #                                  nucleotide.counts[c("chrS2","chrS4","chrS8"),"G"] * Mw_G + 
    #                                  ratio_4sU * nucleotide.counts[c("chrS2","chrS4","chrS8"),"T"] * Mw_4sU + 
    #                                  Mw_triposphate,
    #                              nucleotide.counts[c("chrS5","chrS9","chrS12"),"A"]*Mw_A + nucleotide.counts[c("chrS5","chrS9","chrS12"),"T"]*Mw_U + nucleotide.counts[c("chrS5","chrS9","chrS12"),"C"]*Mw_C + nucleotide.counts[c("chrS5","chrS9","chrS12"),"G"]*Mw_G + Mw_triposphate)
    #for labeled spike-ins
    extended.nucleotide.counts = cbind(nucleotide.counts, nucleotide.counts[,"T"])
    #unlabeled spikeins
    spikein.molecular.weight = c(rowSums(t(t(nucleotide.counts[as.character(seqnames(spikeinCounts))[rowData(spikeinCounts)$labeledSpikein == FALSE], ]) * c(Mw_A, Mw_C, Mw_G, Mw_U))),
       rowSums(t(t(extended.nucleotide.counts[as.character(seqnames(spikeinCounts))[rowData(spikeinCounts)$labeledSpikein == TRUE], ]) * c(Mw_A, Mw_C, Mw_G, (1-ratio_4sU) * Mw_U, ratio_4sU*Mw_4sU)))
       ) + Mw_triposphate
    
    #TODO:
    #calculate spike-in numbers (absolute) (needs amount of spikeins (in nanogram?))
    
    #calcualte spike-in numbers per cell (needs cell number)
}