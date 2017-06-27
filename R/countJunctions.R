#' count Donor, Acceptor and Junction reads
#'
#' @param rows GRanges
#' @param designMatrix dgM 
#' 
#' @return An empty rCubeExperiment container
#' @export
#' @author Leonhard Wachutka
#'
#' @examples
#' 
countJunctions = function(experimentalSetup, filter = NULL, BPPARAM = NULL)
{
	#first count spliced reads
	bamFiles = colData(experimentalSetup)$filename
	
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
		sc[,filename:=filename]
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
	#res = res[,.(count=sum(count)),by=c("filename","seqnames","start","end","strand")]
	
	
}



