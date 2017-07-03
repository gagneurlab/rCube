

estimateRateByFristOrderKineticsSeries = function(featureCounts, rRates, replicate = NULL)
{
	
}

#TODO: move rep into metadata
estimateRateByFristOrderKineticsSeriesCondition = function(featureCounts, replicates, conditions, BPPARAM = BPPARAM,rep = 3 , verbose = FALSE)
{
	#featureCounts = rCubeCounts
	#conditions = 'DMSO.30'
	#replicates = c('1.dedup.bam')
	#BPPARAM = MulticoreParam(workers = 4, progressbar = TRUE)
	fset = subset(featureCounts,,condition==conditions & replicate %in% replicates)
	
	index0 = 1:nrow(fset)
	index1 = lapply(1:rep, function(x)sample(index0))
	batchSize = 100
	
	batches = Reduce(c,lapply(index1, function(x) split(x, ceiling(seq_along(x)/batchSize))))
	batches = batches[1:100]
	
	res = bplapply(batches, callFit, experiment = fset, BPPARAM = BPPARAM, verbose = verbose)
	data = rbindlist(res)
	#Take the median of all refits
	data2 = data[,.(gm = median(gs), gs = median(gs)),by=.(seqnames,start,end,strand,typ)]
	
}
callFit = function(batch, experiment, verbose = FALSE)
{
	if(verbose)
	{
		message('Running on batch ', batch)
	}
	F = colData(experiment)$sizeFactor
	F = (F/F[1])#[-1]
	
	modelTime = as.numeric(colData(experiment)$labelingTime)
	modelTime[colData(experiment)$LT=='T'] = Inf
	ss = experiment[batch]
	
	fit = SE_fit_rates(assay(ss),modelTime,length = rep(1, nrow(ss)),uc = rep(0, nrow(ss)), puc = 0, F = F, gc = 0)
	rr = rowRanges(ss)
	rr$gs = fit$gs
	rr$gm = fit$gm
	return(as.data.table(rr))

}

#the steady states have to be on the end, because of the way we handle gc!!!
#maybe later do a check for this
SE_fit_rates = function(counts,ti, length, uc, puc, params.initial = guess_initial2(counts, ti,...),log.out = FALSE,...)
{
	##Do parametrization
	p = parametrize2(params.initial)
	#Do this so params are in the right order, we cant assume this since there could be user defined params.initial
	params.initial <- list(gs=p$gs,gm=p$gm)
	#ubias = 1-(puc)^uc
	#try with 0 ubias
	ubias = 1
	params.fixed <- list(t=ti,F0 = 1,ga=40,N=nrow(counts), numFinite = sum(is.finite(ti)),l=length,gl = p$gl,gc = p$gc, F = p$F, uc = ubias)
	
	if(log.out != FALSE){print("Start fit...")}
	log_trace = 0
	if(log.out > 0){log_trace = log.out}
	
	fit<-optim(unlist(as.relistable(params.initial)), SElikelihood2, gr=SElikelihood2.grad, params.fixed,params.initial, counts,control=list(maxit=20000,trace=log_trace), method="BFGS")	
	
	if(log.out != FALSE){print("Fit finish...")}
	
	#Undo parametrization of fitted values and collecting them for output
	
	
	params.all = merge_and_repar2(fit$par,params.fixed,params.initial)
	params.all = append(params.all,list(params.initial = params.initial,fit = fit, gen.num = nrow(counts)))
	return (params.all)
}



