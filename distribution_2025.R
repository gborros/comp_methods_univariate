## Testing Stratification Alternative - Dividing up distribution (percentile-based approach)
## Georgi Borros
## 29-06-2025

# Set working directory to path with 'reproducible-code' folder
setwd("C:/Users/01459189/OneDrive/phd/r/sa_comp_methods/reproducible-code")

library(haven)
library(stratification)
library(dplyr)
source("dfstrata.R")
################################################################################
## Part 1: Load Datasets
################################################################################

##### Dataset 1: ECD Census 2021 ######
load("ecd_census.RData")
df <- zap_label(zap_labels(df[order(df$fees_amount),]))
ecd <- as.data.frame(df[, c("fees_amount")])
names(ecd) <- "y"
ecd <- as.data.frame(ecd[ecd$y!=0,])
names(ecd) <- "y"
#ecd_n <- 3711 ## 10% sample size, rounded up
ecd_n <- 20250 ## sample size needed to achieve cv of 0.0037 (comment this out if wanting to see results for 10% sample size and use the line above)


#### Dataset 2: GHS 2024 ####

load("ghs_2024.RData")
df <- zap_label(zap_labels(df[order(df$totmhinc),]))
ghs <- as.data.frame(df[, c("totmhinc")])
names(ghs) <- "y"
#ghs_n <- 2051
ghs_n <- 18000 ## sample size needed to achieve cv of 0.0049

#### Dataset 3: GHS 2023 ####

load("ghs_2023.RData")
df <- zap_label(zap_labels(df[order(df$totmhinc),]))
ghs23 <- as.data.frame(df[, c("totmhinc")])
names(ghs23) <- "y"
#ghs23_n <- 2055
ghs23_n <- 18575 ## sample size needed to achieve cv of 0.0049


#### Dataset 4: TB5 ####
load("tb5.RData")
df <- zap_label(zap_labels(df[df$age<=59,]))
df <- df[order(df$total_elom),]
tb5<- as.data.frame(df[, c("total_elom")])
names(tb5) <- "y"
#tb5_n <- 523
tb5_n <- 640 ## sample size needed for cv of 0.0034

################################################################################
## Part 2: Generate Results for each dataset
################################################################################

###########################ECD Census 5 STRATA#######################################

# data prep
strata_num <- 5
y <- as.vector(ecd$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()

quintiles <- quantile(ecd$y, probs = seq(0, 1, 0.2), na.rm = TRUE)
strata <- cut(ecd$y, breaks = quintiles, include.lowest = TRUE, labels = FALSE)
OSB <- quintiles[2:5]
res <- strata.bh(x=ecd$y, bh=OSB, Ls=5, n=ecd_n, alloc=list(q1=0.5, q2=0, q3=0))

nh <- as.vector(res$nh)
n <- sum(nh)

  # Put results into the dataframe
df <- dfstrata(df=ecd, X=y, h= strata_num, OSB=OSB, nh=nh)
  
  ## Calculate accordingly
count <- aggregate(df$y, list(df$strata), FUN=length)
Nh <- count[,2]
  
  ## Stratum Mean
mean <- aggregate(df$y, list(df$strata), FUN=mean)
meanh <- mean[,2]
  
  ## Mean overall
  mean_init <- 0
  for (l in 1:strata_num) {
    mean_init <- mean_init + Nh[l]*meanh[l]
  }
  mean_est <- mean_init/N 
  
  ## Stratum Variance
  varh <- numeric()
  for (l in 1:strata_num) {
    sigma2 <- 0
    select <- df$y[df$strata == paste0("strata",l)] ## select only obs in the strata
    m <- mean(select) ## mean in selected strata
    length <- length(select) ## obs in strata
    
    for (j in 1:length) {
      sigma2 <- sigma2 + (select[j]- m)^2/length
    }
    varh[l] <- sigma2
  }
  
  var_init <- 0
  ## variance of the estimate from Dalenius & Hodges; Cochran; Er
  for (l in 1:strata_num) {
    
    var_init <- var_init + Nh[l]^2*(varh[l]/nh[l])*(1-(nh[l]/Nh[l]))
  }
  
  var_est <- var_init/N^2
  sd_est <- sqrt(var_est)
  cv_est <- sd_est/mean_est
  
  var_res <- rbind(var_res, var_est)
  cv_res <- rbind(cv_res, cv_est)
  n_res <- rbind(n_res, nh)
  bh_res <- rbind(bh_res, OSB)


precision <- cbind(var_res, cv_res)
ecd_results <- list(precision = precision, 
                    n = n_res,
                    osb = bh_res)

mean(precision[,2])
#save(ecd_results, file="ecd_results_distribution.Rdata")
#save(ecd_results, file="ecd_results_distribution_prop.Rdata")
###########################GHS 2024 5 STRATA##########################################

# data prep
strata_num <- 5
y <- as.vector(ghs$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()

quintiles <- quantile(ghs$y, probs = seq(0, 1, 0.2), na.rm = TRUE)
strata <- cut(ghs$y, breaks = quintiles, include.lowest = TRUE, labels = FALSE)
OSB <- quintiles[2:5]
res <- strata.bh(x=y, bh=OSB, Ls=5, n=ghs_n, alloc=list(q1=0.5, q2=0, q3=0))

nh <- as.vector(res$nh)
n <- sum(nh)
  # Put results into the dataframe
  df <- dfstrata(df=ghs, X=y, h= strata_num, OSB=OSB, nh=nh)
  
  ## Calculate accordingly
  count <- aggregate(df$y, list(df$strata), FUN=length)
  Nh <- count[,2]
  
  ## Stratum Mean
  mean <- aggregate(df$y, list(df$strata), FUN=mean)
  meanh <- mean[,2]
  
  ## Mean overall
  mean_init <- 0
  for (l in 1:strata_num) {
    mean_init <- mean_init + Nh[l]*meanh[l]
  }
  mean_est <- mean_init/N 
  
  ## Stratum Variance
  varh <- numeric()
  for (l in 1:strata_num) {
    sigma2 <- 0
    select <- df$y[df$strata == paste0("strata",l)] ## select only obs in the strata
    m <- mean(select) ## mean in selected strata
    length <- length(select) ## obs in strata
    
    for (j in 1:length) {
      sigma2 <- sigma2 + (select[j]- m)^2/length
    }
    varh[l] <- sigma2
  }
  
  var_init <- 0
  ## variance of the estimate from Dalenius & Hodges; Cochran; Er
  for (l in 1:strata_num) {
    
    var_init <- var_init + Nh[l]^2*(varh[l]/nh[l])*(1-(nh[l]/Nh[l]))
  }
  
  var_est <- var_init/N^2
  sd_est <- sqrt(var_est)
  cv_est <- sd_est/mean_est
  
  var_res <- rbind(var_res, var_est)
  cv_res <- rbind(cv_res, cv_est)
  n_res <- rbind(n_res, nh)
  bh_res <- rbind(bh_res, OSB)

precision <- cbind(var_res, cv_res)
ghs_results <- list(precision = precision, 
                    n = n_res,
                    osb = bh_res)
mean(precision[,2])
#save(ghs_results, file="ghs_results_distribution.Rdata")
#save(ghs_results, file="ghs_results_distribution_prop.Rdata")
###########################TB5 5 STRATA#########################################

# data prep
strata_num <- 5
y <- as.vector(tb5$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()


quintiles <- quantile(tb5$y, probs = seq(0, 1, 0.2), na.rm = TRUE)
strata <- cut(tb5$y, breaks = quintiles, include.lowest = TRUE, labels = FALSE)
OSB <- quintiles[2:5]
res <- strata.bh(x=y, bh=OSB, Ls=5, n=tb5_n, alloc=list(q1=0.5, q2=0, q3=0))

nh <- as.vector(res$nh)
n <- sum(nh)

  # Put results into the dataframe
  df <- dfstrata(df=tb5, X=y, h= strata_num, OSB=OSB, nh=nh)
  
  ## Calculate accordingly
  count <- aggregate(df$y, list(df$strata), FUN=length)
  Nh <- count[,2]
  
  ## Stratum Mean
  mean <- aggregate(df$y, list(df$strata), FUN=mean)
  meanh <- mean[,2]
  
  ## Mean overall
  mean_init <- 0
  for (l in 1:strata_num) {
    mean_init <- mean_init + Nh[l]*meanh[l]
  }
  mean_est <- mean_init/N 
  
  ## Stratum Variance
  varh <- numeric()
  for (l in 1:strata_num) {
    sigma2 <- 0
    select <- df$y[df$strata == paste0("strata",l)] ## select only obs in the strata
    m <- mean(select) ## mean in selected strata
    length <- length(select) ## obs in strata
    
    for (j in 1:length) {
      sigma2 <- sigma2 + (select[j]- m)^2/length
    }
    varh[l] <- sigma2
  }
  
  var_init <- 0
  ## variance of the estimate from Dalenius & Hodges; Cochran; Er
  for (l in 1:strata_num) {
    
    var_init <- var_init + Nh[l]^2*(varh[l]/nh[l])*(1-(nh[l]/Nh[l]))
  }
  
  var_est <- var_init/N^2
  sd_est <- sqrt(var_est)
  cv_est <- sd_est/mean_est
  
  var_res <- rbind(var_res, var_est)
  cv_res <- rbind(cv_res, cv_est)
  n_res <- rbind(n_res, nh)
  bh_res <- rbind(bh_res, OSB)


precision <- cbind(var_res, cv_res)
tb5_results <- list(precision = precision, 
                    n = n_res,
                    osb = bh_res)
mean(precision[,2])
#save(tb5_results, file="tb5_results_distribution.Rdata")
#save(tb5_results, file="tb5_results_distribution_prop.Rdata")
###########################GHS 2023 5 STRATA########################################

# data prep
strata_num <- 5
y <- as.vector(ghs23$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()


quintiles <- quantile(ghs23$y, probs = seq(0, 1, 0.2), na.rm = TRUE)
strata <- cut(ghs23$y, breaks = quintiles, include.lowest = TRUE, labels = FALSE)
OSB <- quintiles[2:5]
res <- strata.bh(x=y, bh=OSB, Ls=5, n=ghs23_n, alloc=list(q1=0.5, q2=0, q3=0))

nh <- as.vector(res$nh)
n <- sum(nh)

  # Put results into the dataframe
  df <- dfstrata(df=ghs23, X=y, h= strata_num, OSB=OSB, nh=nh)
  
  ## Calculate accordingly
  count <- aggregate(df$y, list(df$strata), FUN=length)
  Nh <- count[,2]
  
  ## Stratum Mean
  mean <- aggregate(df$y, list(df$strata), FUN=mean)
  meanh <- mean[,2]
  
  ## Mean overall
  mean_init <- 0
  for (l in 1:strata_num) {
    mean_init <- mean_init + Nh[l]*meanh[l]
  }
  mean_est <- mean_init/N 
  
  ## Stratum Variance
  varh <- numeric()
  for (l in 1:strata_num) {
    sigma2 <- 0
    select <- df$y[df$strata == paste0("strata",l)] ## select only obs in the strata
    m <- mean(select) ## mean in selected strata
    length <- length(select) ## obs in strata
    
    for (j in 1:length) {
      sigma2 <- sigma2 + (select[j]- m)^2/length
    }
    varh[l] <- sigma2
  }
  
  var_init <- 0
  ## variance of the estimate from Dalenius & Hodges; Cochran; Er
  for (l in 1:strata_num) {
    
    var_init <- var_init + Nh[l]^2*(varh[l]/nh[l])*(1-(nh[l]/Nh[l]))
  }
  
  var_est <- var_init/N^2
  sd_est <- sqrt(var_est)
  cv_est <- sd_est/mean_est
  
  var_res <- rbind(var_res, var_est)
  cv_res <- rbind(cv_res, cv_est)
  n_res <- rbind(n_res, nh)
  bh_res <- rbind(bh_res, OSB)


precision <- cbind(var_res, cv_res)
ghs23_results <- list(precision = precision, 
                      n = n_res,
                      osb = bh_res)
mean(precision[,2])
#save(ghs23_results, file="ghs23_results_distribution.Rdata")
#save(ghs23_results, file="ghs23_results_distribution_prop.Rdata")