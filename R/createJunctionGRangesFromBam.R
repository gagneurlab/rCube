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
#' @import BiocGenerics
#' @import data.table
#' @import magrittr
#' 
#'
#' @return Returns a GRanges object with junctions.
#' @author Leonhard Wachutka
#' 
#' @examples
#' # Gencode annotation of MYC gene
#' data(example.exons)
#' constitutive.exons = createConstitutiveFeaturesGRangesFromGRanges(example.exons, 1)
#' @export
createJunctionGRangesFromBam = function(bamFiles, support = 10, ncores=2, BPPARAM = NULL){
	
	extractSplicedReads = function(filename,chr)
	{
		suppressPackageStartupMessages(require(GenomicAlignments,quietly=TRUE))
		suppressPackageStartupMessages(require(data.table,quietly=TRUE))
		suppressPackageStartupMessages(require(magrittr,quietly=TRUE))
		suppressPackageStartupMessages(require(data.table))
		
		#Length is maximum ScanBam can handle
		which = GRanges(seqnames = as.character(chr), IRanges(1, 536870912))
		spliced_reads = readGAlignmentPairs(filename,param = ScanBamParam(which=which, flag=scanBamFlag(isSecondaryAlignment=FALSE)))%>%junctions()%>%reduce()%>%unlist()%>%as.data.table()
		sc = spliced_reads[,.(count = .N),by=c("seqnames","start","end","strand")]
		#message('Loaded: ',filename,'...',chr)
		return(sc)
	}
	suppressPackageStartupMessages(require(data.table,quietly=TRUE))
	suppressPackageStartupMessages(require(GenomicAlignments,quietly=TRUE))
	chrs = unique(unlist(lapply(bamFiles,getChrName)))
	param = data.table(expand.grid(chr = chrs, filename = bamFiles, stringsAsFactors=FALSE))
	
	if(is.null(BPPARAM))
	{
		#BPPARAM = MulticoreParam(workers = ncores)
		BPPARAM= SnowParam(workers = ncores, tasks=nrow(param), type = "SOCK", progressbar = TRUE)
	}
	
	res = rbindlist(bpdtapply(param,extractSplicedReads,BPPARAM = BPPARAM))
	return(GRanges(res[,.(count=sum(count)),by=c("seqnames","start","end","strand")][count>=support,c("seqnames","start","end","strand")]))
}


getChrName = function(file)
{
	names(Rsamtools::scanBamHeader(file)[[1]]$targets)
}
bpdtapply = function(dt,FUN,...,BPREDO = list(), BPPARAM=bpparam())
{
	stopifnot(is.data.table(dt))
	bplapply(split(dt, 1:NROW(dt)),bpdtapply_helper,FUNEXPORT = FUN,BPREDO = BPREDO, BPPARAM = BPPARAM,...)
}

bpdtapply_helper = function(args,FUNEXPORT,...)
{
	do.call(FUNEXPORT,c(args, list(...)))
}

test = function()
{
	require(BiocParallel)
	require(GenomicAlignments)
	require(data.table)
	require(magrittr)
	
	bamFiles = c('/data/ouga03/ag_gagneur/project_local/livia_tt_seq_k562/ProcessedData/Star/HT2000/L1_02_Aligned.sortedByCoord.out.bam',
					'/data/ouga03/ag_gagneur/project_local/livia_tt_seq_k562/ProcessedData/Star/HT2000/L1_05_Aligned.sortedByCoord.out.bam')
	
	mparam = SnowParam(workers = ncores, tasks=nrow(param), type = "SOCK")
	createJunctionGRangesFromBam(bamFiles,BPPARAM = mparam)
}
