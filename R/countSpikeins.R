## author: Leonhard Wachutka

#' @rdname Counting
#' @export
#' @examples
#' data("spikeins")
#' data("spikeinLabeling")
#' folder <- system.file("extdata/Jurkat", package="rCube")
#' bamfiles <- list.files(folder, pattern="*.bam$", full.names=TRUE)
#' spikeinCounts <- setupExperimentSpikeins(rows=spikeins, 
#' files=bamfiles, labelingState=spikeinLabeling)
#' assay(spikeinCounts)
#' spikeinCounts <- countSpikeins(spikeinCounts)
#' assay(spikeinCounts)
countSpikeins <- function(experimentalSetup, scanBamParam=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)), BPPARAM=NULL, verbose=FALSE)
{
    bamFiles <- colData(experimentalSetup)$filename
    region <- rowRanges(experimentalSetup)
    for (fn in bamFiles)
    {
        if(verbose){
            message(date(), ' Counting reads from ', fn)
        }
        
        counts <- .countSpike(fn, region, scanBamParam)
        assays(experimentalSetup[, experimentalSetup[['filename']] == fn])[['counts']] <- as.matrix(counts, ncol=1)
    }
    return(experimentalSetup) 
}


.countSpike <- function(bamFile, region, scanBamParam)
{
    bamWhich(scanBamParam) <- region
    #this will split the Readpairs by CIGAR and merge(union) the resulting reads to avoid double counting of the two ends
    exploded_reads <- GenomicAlignments::readGAlignmentPairs(bamFile, param=scanBamParam)
    exploded_reads <- reduce(grglist(exploded_reads))
    count <- countOverlaps(region, exploded_reads, minoverlap=2)
    return(count)
}


#' @rdname Counting
#' @export
#' @examples
#' data(exampleExons)
#' folder <- system.file("extdata/Jurkat", package="rCube")
#' bamfiles <- list.files(folder, pattern="*.bam$", full.names=TRUE)
#' exonCounts <- setupExperiment(exampleExons, designMatrix=NULL, files=bamfiles)
#' assay(exonCounts)
#' exonCounts <- countFeatures(exonCounts)
#' assay(exonCounts)
countFeatures <- function(experimentalSetup, scanBamParam=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)), BPPARAM=NULL, verbose=FALSE)
{
    countSpikeins(experimentalSetup, scanBamParam=scanBamParam, BPPARAM=BPPARAM, verbose=verbose)
        
}