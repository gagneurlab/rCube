

estimateRateByFristOrderKineticsSeries = function(featureCounts,replicate)
{
	
}

#TODO: move rep into metadata
estimateRateByFristOrderKineticsSeriesCondition = function(featureCounts, replicates, conditions,rep = 3)
{
	#featureCounts = rCubeCounts
	#conditions = 'DMSO.30'
	#replicates = c('1.dedup.bam')
	fset = subset(featureCounts,,condition==conditions & replicate %in% replicates)
	F = colData(fset)$sizeFactor
	F = (F/F[1])#[-1]
	
	modelTime = as.numeric(colData(fset)$labelingTime)
	modelTime[colData(fset)$LT=='T'] = Inf
	
	ss = fset[1:100]
	
	
	SEfit_rates2(assay(ss),modelTime,length = rep(1, nrow(ss)),uc = rep(0, nrow(ss)), puc = 0, F = F, gc = 0)
	
}


#the steady states have to be on the end, because of the way we handle gc!!!
#maybe later do a check for this
#This function fits the model to a bunch of junctions.
SEfit_rates2 = function(counts,ti, length, uc, puc, params.initial = guess_initial2(counts, ti,...),log.out = FALSE,...)
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


