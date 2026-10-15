# Computational Methods for stratified sampling applied to South African data. 

The following files detail how results were generated (all saved within shared folder 'reproducible-code'):

# Prior to running the methods:
We use the following datasets: ghs_2024.RData, ghs_2023.RData, ecd_census.RData, tb5.RData
We use the following author-written functions: elbow_plot.R, dfstrata.R
We use the following as cv input to the n minimising methods: min_cv.RData (generated using script cvmm_packages.R once the cv minimising methods were run)

# Package/method comparison was run as follows:

## stratvns_2025.R: file to run method vns (variable neighbourhood search)
## stratifyR_2025.R: file to run method stratifyR (dynamic programming)
## stratification_2025.R: file to run method rs (random search)
## stratbr_2024.R: file to run method brkga (biased random key GA)
## samplingstrata_2025.R: file to run method ss (grouping genetic algorithm)
## kmeans_2025.R: file to run method kmeans
## ga4stratification_2025.R: file to run GA method
## distribution_2025.R: file to run ptile method

