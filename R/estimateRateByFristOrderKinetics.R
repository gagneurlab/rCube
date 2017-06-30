#' Title
#' 
#' @return rCubeExperiment experiment with rates
#' @export
#' @author Leonhard Wachutka
#'
estimateRateByFristOrderKinetics = function(featureCounts, replicate, method = c('series','single'))
{
	if(method =='series')
	{
		return(estimateRateByFristOrderKineticsSeries(featureCounts,replicate))
	}
}