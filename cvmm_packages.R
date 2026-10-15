## Coefficient of variation-minimising package comparison
## Note: this script is to be run after results have been generated for stratification, stratbr and ga4stratification (cv minimising packages)

# Set working directory to path with 'reproducible-code' folder
setwd("C:/Users/01459189/OneDrive/phd/r/sa_comp_methods/reproducible-code")

library(data.table)
library(haven)
options(scipen=999)

#### Dataset 1: ECD Census ####

### load results obtained
load('ecd_results_stratification.RData')
strat <- ecd_results
load('ecd_results_stratbr.RData')
br <- ecd_results
load('ecd_results_GA4.RData')
ga4 <- ecd_results

### create dataframe that saves the results (variance, cv, runtime by method)
precision <- as.data.frame(rbind(strat$precision, br$precision, ga4$precision))
precision$method <- "rs"
precision$method[11:20] <- "brkga"
precision$method[21:30] <- "ga"
names(precision) <- c("variance", "cv", "runtime", "method")

## Get the minimum cv achieved across all methods and runs
min_ecd <- min(precision$cv)

#### Dataset 2: GHS 2024 ####

### load results obtained
load('ghs_results_stratification.RData')
strat <- ghs_results
load('ghs_results_stratbr.RData')
br <- ghs_results
load('ecd_results_GA4.RData')
ga4 <- ghs_results

### create dataframe that saves the results (variance, cv, runtime by method)
precision <- as.data.frame(rbind(strat$precision, br$precision, ga4$precision))
precision$method <- "rs"
precision$method[11:20] <- "brkga"
precision$method[21:30] <- "ga"
names(precision) <- c("variance", "cv", "runtime", "method")

## Get the minimum cv achieved across all methods and runs
min_ghs <- min(precision$cv)

#### Dataset 3: GHS 2023 ####

### load results obtained
load('ghs23_results_stratification.RData')
strat <- ghs23_results
load('ghs23_results_stratbr.RData')
br <- ghs23_results
load('ghs23_results_GA4.RData')
ga4 <- ghs23_results

### create dataframe that saves the results (variance, cv, runtime by method)
precision <- as.data.frame(rbind(strat$precision, br$precision, ga4$precision))
precision$method <- "rs"
precision$method[11:20] <- "brkga"
precision$method[21:30] <- "ga"
names(precision) <- c("variance", "cv", "runtime", "method")

## Get the minimum cv achieved across all methods and runs
min_ghs23 <- min(precision$cv)

#### Dataset 4: TB5 ####

### load results obtained
load('tb5_results_stratification.RData')
strat <- tb5_results
load('tb5_results_stratbr.RData')
br <- tb5_results
load('tb5_results_GA4.RData')
ga4 <- tb5_results

### create dataframe that saves the results (variance, cv, runtime by method)
precision <- as.data.frame(rbind(strat$precision, br$precision, ga4$precision))
precision$method <- "rs"
precision$method[11:20] <- "brkga"
precision$method[21:30] <- "ga"
names(precision) <- c("variance", "cv", "runtime", "method")

## Get the minimum cv achieved across all methods and runs
min_tb5 <- min(precision$cv)


#### Combine all best results to make one min_cv dataset ####
min_cv <- c(min_ecd, min_tb5, min_ghs, min_ghs23)
names(min_cv) <- c("ecd", "tb5", "ghs", "ghs23")
save(min_cv, file="min_cv.RData")
