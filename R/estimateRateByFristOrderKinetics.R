#' Title
#' 
#' @return rCubeExperiment experiment with rates
#' @export
#' @author Leonhard Wachutka
#'
estimateRateByFristOrderKinetics = function(featureCounts, replicate, method = c('series','single'))
{
	#replicate = list(1,2,c(1,2))
	if(method =='series')
	{
		rRates = createRResultCubeRates(featureCounts,replicate)
		return(estimateRateByFristOrderKineticsSeries(featureCounts,rRates, replicate))
	}
}

createRResultCubeRates = function(featureCounts,replicate)
{
	#designmatrix
	cond = unique(featureCounts$condition)
	dm = expand.grid(condition = cond, rate = as.factor(c('synthesis','decay')), replicate = replicate)
	dm$sample = paste0(dm$condition,'_',dm$rate,'_',dm$replicate)
	
	rates <- matrix(NA, nrow = length(featureCounts), ncol = nrow(dm))
	
	se <- SummarizedExperiment(assays = list("rates"=rates), rowRanges =rowRanges(featureCounts),colData = dm)
	colnames(se) = dm$sample
	se <- new("rCubeRates", se)
	return(se)
	
}