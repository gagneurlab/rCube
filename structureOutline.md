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

# Flow
## Preprocessing
### Create GTF annotation from all extracted junctions (split read)
*In*: list of bams
*Out*: GRanges (GTF)

- (GRanges) createJunctionGRangesFromBam(string[] bamFileVector)
- (GRanges) createJunctionGRangesFromBam(string path)

### Create GTF annotation from constitutive exons/introns (from existing GTF)
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

## Counting
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
Preloaded with default values from typicall TTSeq spiekein set

- (summarizedExperiment) countSpikeins(summarizedExperiment experimentalSetupSpikeins)

## Normalization
*In*: summarized experiment (SE10), full summarized experiment (spikein, SE10s)
*Out*: summarized experiment with additional normalization column (SE20)
- estimate by Carina (spikein)
- estimate by Leo (spikein)
- estimate by Leo (joint model)
	- (summarizedExperiment) calculateNormalization(summarizedExperiment featureCounts, summarizedExperiment spikeinCounts, method=c('spikeinGLM','spikeinMean','jointModel'))

## Dispersion estimation
*In*: full summarized experiment with additional normalization column (SE20)
*Out*: full summarized experiment with additional normalization column and genwise overdispersion Parameter (SE30)
- DESeq
- mean/var by replicates

	- (summarizedExperiment) calculateNormalization(summarizedExperiment featureCounts, method=c('DESeqDispMAP','DESeqDispFit','DESeqDispGeneEst','Replicate'))

## (optional) Model Time
*In/Out*:  SE30
- Add model time column to SE30

	- (summarizedExperiment) calculateModelTime(summarizedExperiment featureCounts, numeric[] time)

## Rate estimation
*In*: SE30
*Out*: Count Rate Matrix

replicate can be *NULL*: In case of 3 replicates this means A, B, C, ABC
Carinas and Leos model get initial parameters from ratio model (only distribution mean and variance)


- ratio (L/T estimation of rates)
	- (summarizedExperiment) calculateRateByRatio(summarizedExperiment featureCounts, factor[] replicate)
	
- Carina
	
	- (summarizedExperiment) calculateRateByCondition(summarizedExperiment featureCounts, factor[] replicate)
	
- Leo
	
	- (summarizedExperiment) calculateRateByLabelingTimeSeries(summarizedExperiment featureCounts, factor[] replicate)

## Postprocessing
*In*: Count Rate Matrix, GTF
*Out*: summarized Rate Matrix
	
	- (summarizedExperiment) mergeRatesByOverlaps(summarizedExperiment featureRates, GRanges topLevelFeatures)
	- (summarizedExperiment) weightedMergeRatesByOverlaps(summarizedExperiment featureRates, GRanges topLevelFeatures, summarizedExperiment featureCounts)

merge different exons/junctions based on findOverlaps





