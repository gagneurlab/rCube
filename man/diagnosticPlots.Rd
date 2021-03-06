% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{plotSpikeinCountsVsSample}
\alias{plotSpikeinCountsVsSample}
\alias{plotResultsReplicates}
\alias{plotFittedCounts}
\title{Diagnostic plots}
\usage{
plotSpikeinCountsVsSample(spikeinCounts)

plotResultsReplicates(featureRates, condition = NULL, replicate = NULL)

plotFittedCounts(featureCounts, featureRates, condition = NULL,
  replicate = NULL)
}
\arguments{
\item{spikeinCounts}{A \code{rCubeExperiment} object containing read 
counts and labeling status for spike-ins, as well as sample information,
see constructor function \code{\link{setupExperimentSpikeins}}.}

\item{featureRates}{A \code{rCubeRates} object containing estimated labeled
and total RNA amounts for all desired replicates/conditions.}

\item{condition}{(optional) Character vector indicating subset of samples
for which plot should be produced}

\item{replicate}{(optional) replicate/combination for which plot should be produced}

\item{featureCounts}{A \code{rCubeExperiment} object containing read counts
and sample information, as well as estimated size factors and 
cross-contamination rates (see \code{\link{estimateSizeFactors}}).}
}
\value{
A \code{trellis} object.
}
\description{
\code{plotSpikeinCountsVsSample} visualizes spike-in read counts
per sample to check if experiment worked and labeled spike-ins are enriched,
but unlabeled spike-ins are depleted in labeled samples, whereas all spike-ins
should be present to a similar extend in total RNA-seq samples.

\code{plotResultsReplicates} plots pairwise comparisons for all
replicates/replicate combinations used for estimating, separated by condition
and rate type from a \code{rCubeExperiment} object.

\code{plotFittedCounts} plots for all investigated conditions and
replicates the observed vs expected read counts in labeled and total RNA-seq
samples.
}
\examples{
data(spikeinCounts)
plotSpikeinCountsVsSample(spikeinCounts)
data(exonRates)
plotResultsReplicates(exonRates)
data(exonCounts)
data(spikeinCounts)
exonCounts <- estimateSizeFactors(featureCounts=exonCounts, 
spikeinCounts=spikeinCounts, method="spikeinGLM")
data(exonRates)
# plot only for ActivatedJurkat condition and one replicate
plotFittedCounts(exonCounts, exonRates, condition="ActivatedJurkat", replicate="1")
# plot for all conditions and estimated replicates/replicate combinations
plotFittedCounts(featureCounts=exonCounts, featureRates=exonRates)
}
\author{
Carina Demel
}
