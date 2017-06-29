#' count Donor, Acceptor and Junction reads
#'
#' @param rows GRanges
#' @param designMatrix dgM 
#' 
#' @import BiocParallel
#' 
#' @return An empty rCubeExperiment container
#' @export
#' @author Leonhard Wachutka
#'
#' @examples
#' 
countJunctions = function(experimentalSetup, filter = NULL, BPPARAM = NULL, ncores = 1, verbose = FALSE)
{
	#first count spliced reads
	bamFiles = colData(experimentalSetup)$filename
	
	suppressPackageStartupMessages(require(data.table,quietly=TRUE))
	suppressPackageStartupMessages(require(GenomicAlignments,quietly=TRUE))
	chrs = unique(seqnames(experimentalSetup))
	param = data.table(expand.grid(chromosome = chrs, bamFile = bamFiles, stringsAsFactors=FALSE))
	
	
	if(is.null(BPPARAM))
	{
		#BPPARAM = MulticoreParam(workers = ncores, progressbar = TRUE)
		BPPARAM = MulticoreParam(workers = ncores, progressbar = TRUE)
		#BPPARAM = SerialParam()
		#BPPARAM = SnowParam(workers = ncores, tasks=nrow(param), type = "SOCK", progressbar = TRUE)
	}
	
	if(verbose){
		message(date(),' Counting J...')
	}
	res = unlist(GRangesList(bpdtapply(param,countSplitReadsPerChromosome,BPPARAM = BPPARAM)))
	resByFilename = split(res,res$filename)
	for (fn in names(resByFilename))
	{
		s = subset(experimentalSetup, typ == 'junction', filename == fn)
		assays(experimentalSetup[rowRanges(experimentalSetup)$typ == 'junction', experimentalSetup[['filename']] == fn])[['counts']] = 
				as.matrix(S4Vectors::merge(rowRanges(s),resByFilename[[fn]], all.x = TRUE)$count,ncol = 1)
	}
	
	
	if(verbose){
		message(date(),' Counting DA...')
	}
	region = subset(rowRanges(experimentalSetup),typ=='donor' | typ=='acceptor')
	res = unlist(GRangesList(bpdtapply(param,countDA,region = region,BPPARAM = BPPARAM)))
	resByFilename = split(res,res$filename)
	require(GenomicAlignments,quietly = TRUE)
	for (fn in names(resByFilename))
	{
		if(verbose){
			message(date(),' Merging file ', fn)
		}
		#merging has to be done in to steps because of unique donor/acceptors
		
		data = subset(resByFilename[[fn]],typ=='donor')
		mCounts = Reduce(function(x,y){S4Vectors::merge(x,y)},split(data,names(data)))
		
		if(verbose){
			message(date(),' Merging donor ', fn)
		}
		assays(experimentalSetup[rowRanges(experimentalSetup)$typ == 'donor', experimentalSetup[['filename']] == fn])[['counts']] = as.matrix(mCounts$count,ncol = 1)
		
		
		data = subset(resByFilename[[fn]],typ=='acceptor')
		if(verbose){
			message(date(),' Merging acceptor ', fn)
		}
		mCounts = Reduce(function(x,y){S4Vectors::merge(x,y)},split(data,names(data)))
		assays(experimentalSetup[rowRanges(experimentalSetup)$typ == 'acceptor', experimentalSetup[['filename']] == fn])[['counts']] = as.matrix(mCounts$count,ncol = 1)
		
	}
	
	return(experimentalSetup)	
}

countDA = function(chromosome, bamFile, region)
{
	suppressPackageStartupMessages(require(GenomicAlignments,quietly = TRUE))
	#message(chromosome)
	which = GRanges(seqnames = chromosome, IRanges(1, 536870912))
	#this will split the Readpairs by CIGAR and merge(union) the resulting reads to avoid double counting of the two ends
	exploded_reads = readGAlignmentPairs(bamFile, param = ScanBamParam(which = which,flag=scanBamFlag(isSecondaryAlignment=FALSE)))
	exploded_reads = reduce(grglist(exploded_reads))
	
	region$count = countOverlaps(region,exploded_reads,minoverlap = 2)
	region$count[region$count==0]=NA
	region$filename = bamFile
	
	return(region)
}

#'
#' counting the split reads per chromosome
#' @noRd
countSplitReadsPerChromosome <- function(chromosome, bamFile){
	suppressPackageStartupMessages(require(GenomicAlignments,quietly=TRUE))
#	message(chromosome)
	# restrict to the chromosome only
	which=GRanges(
			seqnames=chromosome,
			ranges=IRanges(0, 536870912)
	)
	#param <- mergeBamParams(bamParam=scanBamParam(settings), which=which)
	param = ScanBamParam(which=which, flag=scanBamFlag(isSecondaryAlignment=FALSE))
	if(is.null(param)){
		return(GRanges())
	}
	
	# get reads from bam file
	galignment <- readGAlignmentPairs(bamFile, param=param)
	
	# remove the strand information if unstranded data
#	if(!strandSpecific(settings)){
#		strand(galignment) <- "*"
#	}
	
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
			type = 'equal'
	)
	values(junctionsCounts) = cbind(values(junctionsCounts), data.frame(filename = bamFile))
	# sort it and return the GRange object
	return(sort(junctionsCounts))
}



