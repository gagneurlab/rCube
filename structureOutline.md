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
*Out*: GTF
### Create GTF annotation from constitutive exons (from existing GTF)
*In*: GTF
*Out*: GTF

## Create Design Matrix and rowData Matrix
*In*: list of bams or path to scan
*Out*: empty summarized experiment

## Counting
### Features

- Regions (Exons, Introns)
- Junctions (Donor, Acceptor, split reads)
*In*: empty summarized experiment (SE0, GTF)
*Out*: full summarized experiment (SE10)

### Spike-ins
- same as Regions function, additional columns in Design Matrix: spikein length, spikein labelling state (L/T)
- preloaded with default values from typicall TTSeq spiekein set
*In*: empty summarized experiment (SE0, GTFs)
*Out*: full summarized experiment (spikein, SE10s)

## Normalization
*In*: summarized experiment (SE10), full summarized experiment (spikein, SE10s)
*Out*: summarized experiment with additional normalization column (SE20)
- estimate by Carina (spikein)
- estimate by Leo (spikein)
- estimate by Leo (joint model)

## Dispersion estimation
*In*: full summarized experiment with additional normalization column (SE20)
*Out*: full summarized experiment with additional normalization column and genwise overdispersion Parameter (SE30)
- DESeq
- mean/var by replicates

## (optional) Model Time
*In/Out*:  SE30
- Add model time column to SE30

## Rate estimation
*In*: SE30
*Out*: Count Rate Matrix

- ratio (L/T estimation of rates)
- Carina
- Leo

Carinas and Leos model get initial parameters from ratio model (only distribution mean and variance)

## Postprocessing
*In*: Count Rate Matrix, GTF
*Out*: summarized Rate Matrix

merge different exons/junctions based on findOverlaps





