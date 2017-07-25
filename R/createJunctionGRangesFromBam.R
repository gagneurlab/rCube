##TODO CARINA modified import data.table here because of error with shift function in SummarizedExp and data.table

#' Title
#' 
#' @param bamFiles A string Vector containing all bam files
#' @param support min number of supporting reads
#' @param ncores Number of cores available for parallel computation
#' @param BPPARAM see BiocParallel
#' 
#' 
#' @import BiocParallel
#' @import GenomicAlignments
#' @importFrom data.table data.table rbindlist copy is.data.table :=  .N
#' @import magrittr
#' @import GenomicRanges
#' @import IRanges
#' @import Rsamtools
#' 
#'
#' @return Returns a GRanges object with junctions.
#' @author Leonhard Wachutka
#' 
#' @examples
#' # Gencode annotation of MYC gene
#' data(exampleExons)
#' constitutive.exons = createConstitutiveFeaturesGRangesFromGRanges(exampleExons, BPPARAM=NULL, 1)
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
        spliced_reads <- GenomicAlignments::readGAlignmentPairs(bamFile, param=ScanBamParam(which=which, flag=scanBamFlag(isSecondaryAlignment=FALSE))) %>% junctions() %>% reduce() %>% unlist() %>% as.data.table()
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
    res <- rbindlist(bpdtapply(param, extractSplicedReads, BPPARAM=BPPARAM))
    res <- res[, .(count=sum(count)), by=c("seqnames", "start", "end", "strand")][count >= support, c("seqnames", "start", "end", "strand")]
    res <- .jToDA(res)
    res <- unique(res, by=c("seqnames", "start", "end", "strand", "typ"))
    return(sort(GenomicRanges::GRanges(res)))
}


.getChrName <- function(file)
{
    names(Rsamtools::scanBamHeader(file)[[1]]$targets)
}


bpdtapply <- function(dt, FUN, ..., BPREDO=list(), BPPARAM=bpparam())
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
    junctions[, start:=start-1] #TODO Carina asks if this is ok like this, is there not an assignment missing?
    temp <- junctions
    da <- rbind(junctions[, .(start, end=start+1, typ=ifelse(strand == '+', 'donor', 'acceptor')), by=c('seqnames', 'strand')],
                junctions[, .(start=end, end=end+1, typ=ifelse(strand != '+', 'donor', 'acceptor')), by=c('seqnames', 'strand')]
    )
    junctionsIn[, typ := 'junction']
    return(rbind(junctionsIn, da))
}


# .test = function()
# {
#     # TODO Test function needs to be removed here. but pls include example data
#     # require(BiocParallel)
#     # require(GenomicAlignments)
#     # require(data.table)
#     # require(magrittr)
#     
#     bamFiles <- c('/data/ouga03/ag_gagneur/project_local/livia_tt_seq_k562/ProcessedData/Star/HT2000/L1_02_Aligned.sortedByCoord.out.bam',
#                   '/data/ouga03/ag_gagneur/project_local/livia_tt_seq_k562/ProcessedData/Star/HT2000/L1_05_Aligned.sortedByCoord.out.bam')
#     
#     mparam <- SnowParam(workers=ncores, tasks=nrow(param), type="SOCK")
#     createJunctionGRangesFromBam(bamFiles, BPPARAM=mparam)
# }
