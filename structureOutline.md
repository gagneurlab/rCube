# General
- one .R -file per Function; overloaded functions go into same file.
- documentation ROxygen
- rcmd check / bioccheck

 
# Matrix
## Count Design Matrix
| sample name [string] | condition [factor] | L/T [factor] | labeling time [numeric] | replicate [factor] | filename [string]     |
|----------------------|--------------------|--------------|-------------------------|--------------------|-----------------------|
| PlaB.30min_L_10_A    | PlaB.30min         | L            | 10                      | A                  | PlaB.30min_L_10_A.bam |
|                      | ...                | ...          | ...                     | ...                | ...                   |
| PlaB.30min_L_15_B    | PlaB.30min         | L            | 15                      | B                  | PlaB.30min_L_15_B.bam |
|                      | ...                | ...          | ...                     | ...                | ...                   |
| PlaB.30min_T_10_A    | PlaB.30min         | T            | 10                      | A                  | PlaB.30min_T_10_A.bam |
## Count Rate Matrix
| name                    | condition [factor] | rate [factor] | replicate [factor] |
|-------------------------|--------------------|---------------|--------------------|
| PlaB.30min_synthesis_A  | PlaB.30min         | synthesis     | A                  |
| PlaB.30min_decay_A      | PlaB.30min         | decay         | A                  |
| PlaB.30min_synthesis_AB | PlaB.30min         | synthesis     | AB                 |

# Flow Diagram
## Preprocessing
### Create GTF annotation from all extracted junctions (split read) [L]
*In*: list of bams
*Out*: GRanges (GTF)

- (GRanges) createJunctionGRangesFromBam(string[] bamFileVector)
- (GRanges) createJunctionGRangesFromBam(string path)

### Create GTF annotation from constitutive exons/introns (from existing GTF) [C]
*In*: GTF
*Out*: GTF
- (GRanges) createConstituiveFeaturesGRangesFromGRange(GRange gtf)

## Create Design Matrix and rowData Matrix 
*In*: list of bams or path to scan
*Out*: empty summarized experiment

- (summarizedExperiment) setupExperiment(GRanges rows, data.frame designMatrix)
- (summarizedExperiment) setupExperiment(GRanges rows, string[] bamFileVector)
- (summarizedExperiment) setupExperiment(GRanges rows, string path)
- (summarizedExperiment) setupExperimentSpikeins(GRanges rows, string path, numeric[] length, factor[] labelingState)

## Counting [L,Ch]
### Features 
*In*: empty summarized experiment (SE0, GTF)
*Out*: full summarized experiment (SE10)

- Regions (Exons, Introns)
	- (summarizedExperiment) countRegions(summarizedExperiment experimentalSetup)

- Junctions (Donor, Acceptor, split reads)
	- (summarizedExperiment) countJunctions(summarizedExperiment experimentalSetup)


### Spike-ins
*In*: empty summarized experiment (SE0, GTFs)
*Out*: full summarized experiment (spikein, SE10s)

Same as Regions function, additional columns in Design Matrix: spikein length, spikein labelling state (L/T)
Preloaded with default values from typicall TTSeq spiekein set. Data already in repo.

- (summarizedExperiment) countSpikeins(summarizedExperiment experimentalSetupSpikeins)

## Normalization
*In*: summarized experiment (SE10), full summarized experiment (spikein, SE10s)
*Out*: summarized experiment with additional normalization column (SE20)
- estimate by Carina (spikein)
- estimate by Leo (spikein)
- estimate by Leo (joint model)
	- (summarizedExperiment) calculateNormalization(summarizedExperiment featureCounts, summarizedExperiment spikeinCounts, method=c('spikeinGLM','spikeinMean','jointModel'))

(Wrapper by Carina, internal methods by whom it belongs. calculateNormalizationBySpikeinGLM, calculateNormalizationBySpikeinMean, calculateNormalizationByJointModel)
## Dispersion estimation
*In*: full summarized experiment with additional normalization column (SE20)
*Out*: full summarized experiment with additional normalization column and genwise overdispersion Parameter (SE30)
Adds two columns (dispersionLabel, dispersionTotal)
- DESeq
- mean/var by replicates

	- (summarizedExperiment) calculateNormalization(summarizedExperiment featureCounts, method=c('DESeqDispMAP','DESeqDispFit','DESeqDispGeneEst','Replicate'))
(wrapper (df with 2 cols) calculateNormalizationBy...Method...)
## (optional) Model Time [L]
*In/Out*:  SE30
- Add model time column to SE30

	- (summarizedExperiment) calculateModelTime(summarizedExperiment featureCounts, numeric[] time)

## Rate estimation
*In*: SE30
*Out*: Count Rate Matrix

replicate can be *NULL*: In case of 3 replicates this means A, B, C, ABC
Carinas and Leos model get initial parameters from ratio model (only distribution mean and variance)


- ratio (L/T estimation of rates) [L]
	- (summarizedExperiment) calculateRateByRatio(summarizedExperiment featureCounts, factor[] replicate)
	
- Carina [C]
	- (summarizedExperiment) calculateRateByCondition(summarizedExperiment featureCounts, factor[] replicate)
	
- Leo [L]
	- (summarizedExperiment) calculateRateByLabelingTimeSeries(summarizedExperiment featureCounts, factor[] replicate)

## Postprocessing [C]
*In*: Count Rate Matrix, GTF
*Out*: summarized Rate Matrix

	- (summarizedExperiment) mergeRatesByOverlaps(summarizedExperiment featureRates, GRanges topLevelFeatures)
	- (summarizedExperiment) weightedMergeRatesByOverlaps(summarizedExperiment featureRates, GRanges topLevelFeatures, summarizedExperiment featureCounts)

merge different exons/junctions based on findOverlaps



# Debug
- Plot: sample vs. spikein counts
- Plot: RepA vs RepB
- Plot: model time vs normalization factor
- [Print: rate estimation by ratio summary]

