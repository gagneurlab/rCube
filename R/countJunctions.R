## author: Leonhard Wachutka

#' @title Counting reads for features
#' 
#' @description \code{countSpikeins} provides a framework to count reads for 
#' artificial spikeins, \code{countJunctions} to count reads for junctions.
#'
#' @param experimentalSetup An empty \code{rCubeExperiment} object for 
#' spikeins/junctions, see \code{\link{setupExperimentSpikeins}} and 
#' \code{\link{setupExperiment}}.
#' @param scanBamParam Parameter to specify which fields are imported from a BAM
#' file.
#' @param BPPARAM An instance of a \code{BiocParallelParam} class, e.g., 
#' \code{\link{MulticoreParam}}, \code{\link{SnowParam}}, 
#' \code{\link{DoparParam}}.
#' @param ncores Number of cores to be used for parallel computation.
#' @param verbose If true, status messages are printed.
#' 
#' @seealso \code{\link[Rsamtools]{ScanBamParam}}
#' @seealso \code{\link{BiocParallelParam}}
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import Rsamtools
#' @import S4Vectors
#' 
#' @return An updated \code{rCubeExperiment} container, with spike-in read 
#' counts.
#' @export
#' @author Leonhard Wachutka
#' @rdname Counting
#'
#' @examples
#' data(junctions)
#' bamfiles <- list.files(system.file("extdata/TimeSeriesExample/", package='rCube'), 
#' pattern=".bam$", full.names=TRUE)
#' # first setup experiment
#' junctionCounts <- setupExperiment(junctions, files=bamfiles)
#' # then count reads on junctions
#' junctionCounts <- countJunctions(junctionCounts, 
#' BPPARAM=BiocParallel::MulticoreParam(1))
countJunctions <- function(experimentalSetup, scanBamParam=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)), BPPARAM=NULL, ncores=1, verbose=FALSE)
{
    #first count spliced reads
    bamFiles <- colData(experimentalSetup)$filename
    
    chrs <- unique(seqnames(experimentalSetup))
    param <- data.table(expand.grid(chromosome=chrs, bamFile=bamFiles,
                                    stringsAsFactors=FALSE))
    
    if(is.null(BPPARAM)){
        #BPPARAM <- MulticoreParam(workers=ncores, progressbar=TRUE)
        BPPARAM <- MulticoreParam(workers=ncores, progressbar=TRUE, log=FALSE)
        #BPPARAM <- SerialParam()
        #BPPARAM <- SnowParam(workers=ncores, tasks=nrow(param), type="SOCK", progressbar=TRUE)
    }
    
    if(verbose){
        message(date(),' Counting J...')
    }
    bptasks(BPPARAM) <- nrow(param)
    res <- unlist(GRangesList(.bpdtapply(param, .countSplitReadsPerChromosome,
                                        scanBamParam=scanBamParam, 
                                        BPPARAM=BPPARAM)))
    resByFilename <- split(res, res$filename)
    for (fn in names(resByFilename))
    {
        s <- subset(experimentalSetup, typ == 'junction', filename == fn)
        assays(experimentalSetup[rowRanges(experimentalSetup)$typ == 'junction',
                                experimentalSetup[['filename']] == fn])[['counts']] <- 
            as.matrix(merge(rowRanges(s), resByFilename[[fn]], all.x=TRUE)$count, ncol=1)
    }
    
    if(verbose){
        message(date(),' Counting DA...')
    }
    
    buildIndex <- function(maxIndex, stepSize=5000)
    {
        end <- (1:ceiling(maxIndex/stepSize)) * stepSize
        start <- end - stepSize + 1
        end[length(end)] <- maxIndex
        data.frame(start, end)
    }
    
    region <- subset(rowRanges(experimentalSetup), 
                    typ == 'donor' | typ == 'acceptor')
    param <- data.table(buildIndex(length(region), stepSize=500))
    #param = param[1:10]
    for (fn in bamFiles)
    {
        if(verbose){
            message(date(),' Counting da for ', fn)
        }
        bptasks(BPPARAM) <- nrow(param)
        res <- .bpdtapply(param, .countDA, bamFile=fn, region=region, 
                        scanBamParam=scanBamParam, BPPARAM=BPPARAM)
        if(verbose){
            message(date(),' Reduce da for ', fn)
        }
        counts  <- Reduce(pmax.int, res)
        if(verbose){
            message(date(),' Assign da for ', fn)
        }
        assays(experimentalSetup[rowRanges(experimentalSetup)$typ == 'donor' | 
                                rowRanges(experimentalSetup)$typ == 'acceptor', 
                                experimentalSetup[['filename']] == fn])[['counts']] <- as.matrix(counts, ncol=1)
    }
    
    # Zero out NA lines.
    assay(experimentalSetup)[is.na(assay(experimentalSetup))] <- 0
    return(experimentalSetup)
}


.countDA <- function(start, end, bamFile, region, scanBamParam)
{
    region2 <- region
    strand(region2) <- '*'
    Rsamtools::bamWhich(scanBamParam) <- range(sort(region2)[start:end])
    #this will split the Readpairs by CIGAR and merge(union) the resulting reads to avoid double counting of the two ends
    exploded_reads <- readGAlignmentPairs(bamFile, param=scanBamParam)
    exploded_reads <- reduce(grglist(exploded_reads))
    count <- countOverlaps(region, exploded_reads, minoverlap=2)
    return(count)
}


#' counting the split reads per chromosome
#' @noRd
.countSplitReadsPerChromosome <- function(chromosome, bamFile, scanBamParam){
    # restrict to the chromosome only
    
    Rsamtools::bamWhich(scanBamParam) <- GRanges(seqnames=chromosome, ranges=IRanges(0, 536870912))
    
    # get reads from bam file
    galignment <- readGAlignmentPairs(bamFile, param=scanBamParam)
    
    # remove the strand information if unstranded data
    #if(!strandSpecific(settings)){
    #   strand(galignment) <- "*"
    #}
    
    # dont count if there is nothing to count
    if(length(galignment) == 0){
        return(GRanges())
    }
    
    # get the junction positions and their counts
    junctions <- unlist(junctions(galignment))
    junctionsCounts <- unique(junctions)
    
    # dont count anything if there is nothing to count
    if(length(junctionsCounts) == 0){
        return(junctionsCounts)
    }
    
    # count the data
    mcols(junctionsCounts)$count <- countOverlaps(
        junctionsCounts,
        junctions,
        type='equal'
    )
    values(junctionsCounts) <- cbind(values(junctionsCounts), 
                                    data.frame(filename=bamFile))
    # sort it and return the GRange object
    return(sort(junctionsCounts))
}
