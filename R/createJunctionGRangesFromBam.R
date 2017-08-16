## author: Leonhard Wachutka
 
#' Create de novo junction Annotation from bam files
#' 
#' @description Reads coming from sequenced fragments that cannot be 
#' mapped entirely to a genomic location, can be split, so that they e.g. map 
#' two exons flanking an intron. The cigar string of these reads usually contains
#' a 'N' character, where the read cannot be aligned to the genome. This information
#' can be used to extract possible splice sites.
#' 
#' @param bamFiles A string vector containing all bam files.
#' @param support Minimal number of supporting reads.
#' @param ncores Number of cores available for parallel computation.
#' @param BPPARAM An instance of a \code{BiocParallelParam} class, e.g., 
#' \code{\link{MulticoreParam}}, \code{\link{SnowParam}}, \code{\link{DoparParam}}.
#' 
#' @seealso \code{\link{BiocParallelParam}}
#' @import BiocParallel
#' @import GenomicAlignments
#' @importFrom data.table data.table rbindlist copy is.data.table :=  .N
#' @import magrittr
#' @import GenomicRanges
#' @import IRanges
#' @import Rsamtools
#' 
#' @return Returns a GRanges object with junctions.
#' @author Leonhard Wachutka
#' 
#' @examples
#' bamfiles <- list.files(system.file("extdata/TimeSeriesExample/", package='rCube'), 
#' pattern=".bam$", full.names=TRUE)
#' # do not run because of time reasons
#' # junctions <- createJunctionGRangesFromBam(bamfiles, support=5, ncores=1, 
#' # BPPARAM=BiocParallel::MulticoreParam(1))
#' @export
createJunctionGRangesFromBam <- function(bamFiles, 
                                         support=10,
                                         ncores=2,
                                         BPPARAM=NULL)
{
    
    extractSplicedReads <- function(bamFile, chromosome)
    {
        #Length is maximum ScanBam can handle
        which <- GRanges(seqnames=as.character(chromosome), IRanges(1, 536870912))
        spliced_reads <- GenomicAlignments::readGAlignmentPairs(bamFile, param=Rsamtools::ScanBamParam(which=which, flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE))) %>% junctions() %>% reduce() %>% unlist() %>% as.data.table()
        sc <- spliced_reads[, .(count= .N), by=c("seqnames", "start", "end", "strand")]
        #message('Loaded: ', bamFile, '...', chromosome)
        return(sc)
    }
    chrs <- unique(unlist(lapply(bamFiles, .getChrName)))
    param <- data.table::data.table(expand.grid(chromosome=chrs, 
                                                bamFile=bamFiles,
                                                stringsAsFactors=FALSE))
    
    if(is.null(BPPARAM))
    {
        BPPARAM <- BiocParallel::MulticoreParam(workers=ncores, progressbar=TRUE)
        #BPPARAM <- SerialParam()
        #BPPARAM <- SnowParam(workers=ncores, tasks=nrow(param), type="SOCK", progressbar=TRUE)
    }
    bptasks(BPPARAM) <- nrow(param)
    res <- rbindlist(.bpdtapply(param, extractSplicedReads, BPPARAM=BPPARAM))
    res <- res[, .(count=sum(count)), by=c("seqnames", "start", "end", "strand")][count >= support, c("seqnames", "start", "end", "strand")]
    res <- .jToDA(res)
    res <- unique(res, by=c("seqnames", "start", "end", "strand", "typ"))
    return(sort(GenomicRanges::GRanges(res)))
}


.getChrName <- function(file)
{
    names(Rsamtools::scanBamHeader(file)[[1]]$targets)
}


.bpdtapply <- function(dt, FUN, ..., BPREDO=list(), BPPARAM=bpparam())
{
    stopifnot(is.data.table(dt))
    bplapply(split(dt, 1:NROW(dt)), .bpdtapply_helper, FUNEXPORT=FUN, BPREDO=BPREDO, BPPARAM=BPPARAM, ...)
}


.bpdtapply_helper <- function(args, FUNEXPORT, ...)
{
    do.call(FUNEXPORT, c(args, list(...)))
}


.jToDA <- function(junctionsIn)
{
    junctions <- copy(junctionsIn)
    junctions[, start:=start-1]
    temp <- junctions
    da <- rbind(junctions[, .(start, end=start+1, typ=ifelse(strand == '+', 'donor', 'acceptor')), by=c('seqnames', 'strand')],
                junctions[, .(start=end, end=end+1, typ=ifelse(strand != '+', 'donor', 'acceptor')), by=c('seqnames', 'strand')]
    )
    junctionsIn[, typ := 'junction']
    return(rbind(junctionsIn, da))
}