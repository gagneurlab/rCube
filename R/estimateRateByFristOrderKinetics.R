#' Title
#' 
#' @return rCubeExperiment experiment with rates
#' @export
#' @author Leonhard Wachutka
#'
estimateRateByFristOrderKinetics = function(featureCounts, replicate, method = c('series','single'), BPPARAM = NULL)
{
	#replicate = c(1,2,'1:2')
	if(method =='series')
	{
		rRates = createRResultCubeRates(featureCounts,replicate)
		return(estimateRateByFristOrderKineticsSeries(featureCounts,rRates, BPPARAM = BPPARAM))
	}else{
	    return(estimateRateByFristOrderKineticsSingle(featureCounts,replicate))
	}
}

createRResultCubeRates = function(featureCounts,replicate)
{
	#designmatrix
	cond = unique(featureCounts$condition)
	dm = expand.grid(condition = cond, rate = as.factor(c('synthesis','degradation')), replicate = as.character(replicate))
	dm$sample = paste0(dm$condition,'_',dm$rate,'_',dm$replicate)
	
	rates <- matrix(NA, nrow = length(featureCounts), ncol = nrow(dm))
	
	se <- SummarizedExperiment(assays = list("rates"=rates), rowRanges =rowRanges(featureCounts),colData = dm)
	colnames(se) = dm$sample
	se <- new("rCubeRates", se)
	return(se)
	
}