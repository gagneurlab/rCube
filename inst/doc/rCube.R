## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE, dev="png", fig.show="hide",
               fig.width=4, fig.height=4.5,
               message=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"---------------------------------------
BiocStyle::latex()

## ----options,results="hide",echo=FALSE----------------------------------------
options(digits=3, width=80, prompt=" ", continue=" ")

## ----LoadingLibrary, echo=TRUE, results="hide"--------------------------------
library("rCube")

## ----importGTF, echo=TRUE, eval=FALSE-----------------------------------------
#  library(rtracklayer)
#  granges <-import(gtffile)

## ----NumberofIterations, eval=FALSE-------------------------------------------
#  elementMetadata(featureCounts)$numberOfIterations <- 10

## ----replicate, eval=FALSE----------------------------------------------------
#  estimateRateByFirstOrderKinetics(featureCounts,
#                                   replicate=c(1, 2, "1:2"),
#                                   method=c("single", "series"),
#                                   BPPARAM=NULL)

## ----constitutiveExons, echo=TRUE---------------------------------------------
data("exampleExons")
exampleExons
constitutiveExons <- createConstitutiveFeaturesGRangesFromGRanges(exampleExons,
                                                                  BPPARAM=NULL,
                                                                  ncores=1)
constitutiveExons

## ----plotConstExons, echo=FALSE, eval=TRUE, fig.width=7.5, fig.height=2.5, fig.show='hide', fig.cap="Illustration of two transcript isoforms for the FOS gene and the resulting constitutive exons"----
library(ggbio)
data("exampleExons")
elementMetadata(exampleExons) = data.frame("group"=elementMetadata(exampleExons)[,"transcript_id"])
elementMetadata(constitutiveExons) = data.frame("group"="constitutive")
gra2 = c(exampleExons, constitutiveExons)
levels(gra2$group) = c("isoform 1", "isoform 2", "constitutive")
gra2$group = relevel(gra2$group, ref="constitutive")
autoplot(gra2, aes(fill = group, group = group), geom = "alignment", group.selfish = TRUE, xlab="chromosome 14")

## ----ExperimentalDesignTcellActivation, echo=TRUE, eval=TRUE------------------
folder <- system.file("extdata/TcellActivation", package='rCube')
expDesign <- read.delim(file.path(folder, "experimentalDesign.txt"))
expDesign

## ----ExperimentalDesignFromFile, echo=TRUE------------------------------------
exonCounts <- setupExperiment(constitutiveExons, designMatrix=expDesign, files=NULL)
class(exonCounts)

## ----bamfiles, echo=TRUE------------------------------------------------------
bamfiles <- list.files(folder, pattern="*.bam$", full.names=TRUE)
basename(bamfiles)

exonCounts <- setupExperiment(constitutiveExons, designMatrix=NULL, files=bamfiles)

## ----expDesignAccessors, echo=TRUE, eval=FALSE--------------------------------
#  # feature information
#  rowRanges(exonCounts)
#  
#  # sample information
#  colData(exonCounts)
#  
#  # read counts
#  assay(exonCounts)

## ----Counting, echo=TRUE, eval=TRUE-------------------------------------------
assay(exonCounts)

exonCounts <- countFeatures(exonCounts)

assay(exonCounts)

## ----spikeins, echo=TRUE------------------------------------------------------
data("spikeins")
data("spikeinLabeling")
spikeinLengths <- width(spikeins)

## ----spikeindesign, echo=TRUE-------------------------------------------------
spikeinCounts <- setupExperimentSpikeins(rows=spikeins,
                                         designMatrix=expDesign,
                                         length=spikeinLengths,
                                         labelingState=spikeinLabeling)

## ----spikeincountsfilename, echo=FALSE----------------------------------------
colData(spikeinCounts)$filename <- bamfiles

## ----spikeincounts, echo=TRUE-------------------------------------------------
spikeinCounts <- countSpikeins(spikeinCounts)
assay(spikeinCounts)

## ----SpikeinDiagnosticPlot, fig.show='hide', fig.width=10, fig.height=4-------
plotSpikeinCountsVsSample(spikeinCounts)

## ----SpikeinNormalization, echo=TRUE, eval=TRUE-------------------------------
exonCounts <- estimateSizeFactors(exonCounts, spikeinCounts, method="spikeinGLM")
colnames(colData(exonCounts))
exonCounts$sizeFactor
exonCounts$crossContamination

## ----SpikeinNormalizationMetadata, echo=TRUE, eval=FALSE----------------------
#  metadata(exonCounts)

## ----DESeqDispersion, echo=TRUE, eval=TRUE------------------------------------
exonCounts <- estimateSizeDispersions(exonCounts, method='DESeqDispGeneEst')
rowRanges(exonCounts)

## ----NumberOfIterations, echo=TRUE, eval=TRUE---------------------------------
elementMetadata(exonCounts)$numberOfInterations <- 7

## ----SynDecEstimation, echo=TRUE, eval=TRUE-----------------------------------
rates <- estimateRateByFirstOrderKinetics(exonCounts,
                                          replicate=c(1, 2, "1:2"),
                                          method='single',
                                          BPPARAM=BiocParallel::MulticoreParam(1))
rates

## ----ObservedVsFittedCounts, fig.show='hide', fig.width=10, fig.height=4------
plotFittedCounts(exonCounts, rates)

## ----ReplicateResults, fig.show='hide', fig.width=10, fig.height=4------------
plotResultsReplicates(rates)

## ----SummarizeRates, echo=TRUE, eval=TRUE-------------------------------------
topLevelFeature <- GRanges(seqnames="chr14", 
                           ranges=IRanges(start=75278774, end=75281636), 
                           strand="+")
topLevelFeaturesRates <- summarizeRates(rates, topLevelFeature, by='mean')

## ----adjustSynthesis, echo=TRUE, eval=TRUE------------------------------------
data("spikeinGenome")
adjustedSynthesis <- adjustSynthesisRate(spikeinCounts=spikeinCounts, 
                                         spikeinGenome=spikeinGenome, 
                                         featureRates=rates, 
                                         spikeinAmount=25*0.05, 
                                         cellNumber=5*10^7)
range(adjustedSynthesis)

## ----junctionAnno, echo=TRUE, eval=FALSE--------------------------------------
#  junctions <- createJunctionGRangesFromBam(bamfiles, support=5, ncores=1,
#                                            BPPARAM=BiocParallel::MulticoreParam(1))

## ----loadJunctions, echo=TRUE, eval=TRUE--------------------------------------
data(junctions)

## ----showJunctions, echo=TRUE, eval=TRUE--------------------------------------
head(junctions)

## ----ExperimentalDesignTimeSeriesExample, echo=TRUE, eval=TRUE----------------
folder <- system.file("extdata/TimeSeriesExample", package='rCube')
expDesign <- read.delim(file.path(folder, "experimentalDesign.txt"))
head(expDesign)

## ----ExperimentalDesignFromFileTimeSeriesExample, echo=TRUE, eval=TRUE--------
junctionCounts <- setupExperiment(junctions, designMatrix=expDesign)
class(junctionCounts)

## ----CountingTimeSeriesExample, echo=TRUE, eval=FALSE-------------------------
#  junctionCounts <- countJunctions(junctionCounts, BPPARAM=BiocParallel::MulticoreParam(1))

## ----assignCounts, echo=TRUE, eval=TRUE---------------------------------------
data(countMatrixJunctions)
assay(junctionCounts) <- countMatrixJunctions

## ----spikeins2, echo=TRUE-----------------------------------------------------
data("spikeins")
data("spikeinLabeling")
spikeinLengths <- width(spikeins)

## ----spikeindesign2, echo=TRUE------------------------------------------------
spikeinCounts <- setupExperimentSpikeins(rows=spikeins,
                                         designMatrix=expDesign,
                                         length=spikeinLengths,
                                         labelingState=spikeinLabeling)

## ----spikeincountsTimeseries, echo=TRUE, eval=FALSE---------------------------
#  spikeinCounts <- countSpikeins(spikeinCounts)

## ----spikeincountsTimeseries2, echo=TRUE, eval=TRUE---------------------------
data("countMatrixSpikeins")
assay(spikeinCounts) <- countMatrixSpikeins

## ----SpikeinNormalization2, echo=TRUE, eval=TRUE------------------------------
junctionCounts <- estimateSizeFactors(junctionCounts, spikeinCounts, method="spikeinMean")
colnames(colData(junctionCounts))
junctionCounts$sizeFactor

## ----Dispersion, echo=TRUE, eval=TRUE-----------------------------------------
junctionCounts <- estimateSizeDispersions(junctionCounts, method='Replicate')
rowRanges(junctionCounts)

## ----Filtering, echo=TRUE, eval=TRUE------------------------------------------
junctionCounts <- junctionCounts[rowSums(assay(junctionCounts)) > 100]

## ----Fitting, echo=TRUE, eval=TRUE--------------------------------------------
junctionRates <- estimateRateByFirstOrderKinetics(junctionCounts,
                                                 replicate=c(1,2,'1:2'),
                                                 method='series')

## ----SummarizeRates2, echo=TRUE, eval=TRUE------------------------------------
topLevelFeaturesJunctions <- GRanges(seqnames="chr3",
                           ranges=IRanges(start=169769500, end=169789800),
                           strand="+")
topLevelFeaturesRates <- summarizeRates(junctionRates, topLevelFeaturesJunctions, by='median')

assay(topLevelFeaturesRates)

## ----halflive, echo=TRUE, eval=TRUE-------------------------------------------
hl <- log(2)/assay(topLevelFeaturesRates)[, topLevelFeaturesRates$rate == "degradation"]
hl

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

