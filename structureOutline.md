# Matrix
## Design Matrix
| sample name [string] | condition [factor] | L/T [factor] | labeling time [numeric] | replicate [factor] | filename [string]     |
|----------------------|--------------------|--------------|-------------------------|--------------------|-----------------------|
| PlaB.30min_L_10_A    | PlaB.30min         | L            | 10                      | A                  | PlaB.30min_L_10_A.bam |
|                      | ...                | ...          | ...                     | ...                | ...                   |
| PlaB.30min_L_15_B    | PlaB.30min         | L            | 15                      | B                  | PlaB.30min_L_15_B.bam |
|                      | ...                | ...          | ...                     | ...                | ...                   |
| PlaB.30min_T_10_A    | PlaB.30min         | T            | 10                      | A                  | PlaB.30min_T_10_A.bam |

# Flow
## Preprocessing
### Create GTF annotation from all extracted junctions (split read)
### Create GTF annotation from constitutive exons (from existing GTF)
## Counting
### Features
- Regions (Exons, Introns)
- Junctions (Donor, Acceptor, split reads)

### Spike-ins

## Normalization

## Dispersion estimation
## (optional) Model Time
## Rate estimation
## Postprocessing







