#' Title
#' 
#' @return rCubeExperiment Object with updated metadata information
#' @export
#' @author Leonhard Wachutka
#'
estimateSizeDispersions = function(experiment, method = c('DESeqDispMAP','DESeqDispFit','DESeqDispGeneEst','Replicate'))
{
	if(method=='Replicate')
	{
		return(estimateSizeDispersions_replicate(experiment))
	}
}

estimateSizeDispersions_replicate = function(experiment)
{
	rowRanges(experiment)$dispersion = 40
	return(experiment)
}