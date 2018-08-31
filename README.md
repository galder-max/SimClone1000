# Simclone grid design

## Installation
There is no installation required.
### Dependencies
R libraries _parallel_ and _simclone_ [github.com/sdentro/simclone](github.com/sdentro/simclone)  
A slurm-based cluster to submit the runs as jobs in _submitSimulateAll.R_

## Run pipelines in R
### Step1
run _deriveGrid.R_ to generate the grid design:
> Rscript derivegrid.R

### Step2
run _submitSimulateAll.R_ to simulate all tumours from the grid design (require simclone library):
> Rscript submitSimulateAll.R

This will submit one job per simulated sample on a slurm-based cluster.

## Output
Each sample run will write the outputs generate by the _simclone_ main function [github.com/sdentro/simclone](github.com/sdentro/simclone)  
