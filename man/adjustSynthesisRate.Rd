% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpikeinAdjustment.R
\name{adjustSynthesisRate}
\alias{adjustSynthesisRate}
\title{Adjust synthesis rate to transcripts per cell per minute}
\usage{
adjustSynthesisRate(spikeinCounts, spikeinGenome, featureRates,
  ratio_4sU = 0.1, spikeinAmount = 1, cellNumber = 10^7)
}
\arguments{
\item{spikeinCounts}{\code{rCubeExperiment} object with spikeins read counts}

\item{spikeinGenome}{A \code{DNAStringSet} containing the spike-in sequences.
Names must be identical to the \code{seqnames} in the \code{spikeinCounts} object.}

\item{featureRates}{\code{rCubeRates} object with resulting synthesis rates 
that should be adjusted.}

\item{ratio_4sU}{Percentage of 4sU vs U volume for reverse transcription of 
spike-ins.}

\item{spikeinAmount}{amount of spike-ins added to each sample in ng. Either 
a single value, if same amounts were used for all spike-ins, or a vector.
Alternatively, a column in the \code{rowData} of the \code{spikeinCounts} can
be set, named 'spikeinAmount'.}

\item{cellNumber}{Number of cells used for the experiment. Either a single 
value if same cell numbers were used for all samples, or a vector.
Alternatively, a column in the \code{colData} of the \code{spikeinCounts} can
be set, named 'cellNumber'.
Note: if you want to use different spike-ins amounts for different samples
with different cell numbers,
please split the data according to the samples and run it separately.}
}
\value{
Returns a matrix of synthesis rates, normalized to 
transcripts/cell/minute.
}
\description{
Based on the knowledge how many spike-ins were originally added 
to a defined
number of cells, estimated synthesis rates can be adjusted to retrieve
a rate for transcripts per cell and minute.
}
\examples{
data("spikeinGenome")
data("spikeinCounts")
# due to downsampling of read counts to 5\% of originial counts, we need to adjust here:
rowRanges(spikeinCounts)$spikeinAmount <- 25*0.05
colData(spikeinCounts)$cellNumber <- 5*10^7
data("exonRates")
adjustedSynthesis <- adjustSynthesisRate(spikeinCounts, spikeinGenome, exonRates)
adjustedSynthesis
}
\author{
Bj\"orn Schwalb, Carina Demel
}
