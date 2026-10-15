## Variable neighbourhood search algorithm package comparison testing
## Author: Georgi Borros
## Date: 12/06/2025

### Note to user: must install stratvns_1.1.tar (available in reproducible-code folder) prior to running this file

# Set working directory to path with 'reproducible-code' folder
setwd("C:/Users/01459189/OneDrive/phd/r/sa_comp_methods/reproducible-code")

library(sampling)
library(data.table)
library(stratification)
library(factoextra)
library(NbClust)
library(haven)
options(scipen=999)
library(stratvns)
library(kableExtra)

# Source author-written functions
source("dfstrata.R")
source("elbow_plot.R")
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
ecd_n <- 3711 ## 10% sample size, rounded up

#### Dataset 2: GHS 2024 ####

load("ghs_2024.RData")
df <- zap_label(zap_labels(df[order(df$totmhinc),]))
ghs <- as.data.frame(df[, c("totmhinc")])
names(ghs) <- "y"
ghs_n <- 2051

#### Dataset 3: GHS 2023 ####

load("ghs_2023.RData")
df <- zap_label(zap_labels(df[order(df$totmhinc),]))
ghs23 <- as.data.frame(df[, c("totmhinc")])
names(ghs23) <- "y"
ghs23_n <- 2055

#### Dataset 4: TB5 ####
load("tb5.RData")
df <- zap_label(zap_labels(df[df$age<=59,]))
df <- df[order(df$total_elom),]
tb5<- as.data.frame(df[, c("total_elom")])
names(tb5) <- "y"
tb5_n <- 523


#### Load minimum cvs achieved by cv minimising packages:
## Min
load('min_cv.RData')

################################################################################
## Part 2: Generate Results for each dataset
################################################################################

###########################ECD Census 5 STRATA##################################

# data prep
strata_num <- 5
y <- as.vector(ecd$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
min_ecd <- as.numeric(min_cv[1])

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  res <- STRATVNS(ecd$y, L=strata_num, cvt=0.0036571, imax=5, maxstart=2, tmax=5, kmax = 3, nIterWithNoImpMax=2, nsols = 10, parallelize=TRUE, cputime=1000)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  
  nh <- as.vector(res$nh)
  n <- sum(nh)
  OSB <- as.vector(res$bk)
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
  time_res <- rbind(time_res, time)
  
}


precision <- cbind(var_res, cv_res, time_res)
ecd_results <- list(precision = precision, 
                    n = n_res,
                    osb = bh_res)

save(ecd_results, file="ecd_results_stratvns.Rdata")

###########################GHS 2024 5 STRATA##########################################

# data prep
strata_num <- 5
y <- as.vector(ghs$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
min_ghs <- as.numeric(min_cv[3])

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  res <- STRATVNS(ghs$y, L=strata_num, cvt=0.0048608, imax=5, maxstart=2, tmax=5, kmax = 3, nIterWithNoImpMax=2, nsols = 10, parallelize=TRUE, cputime=1000)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  
  nh <- as.vector(res$nh)
  n <- sum(nh)
  OSB <- as.vector(res$bk)
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
  time_res <- rbind(time_res, time)
  
}


precision <- cbind(var_res, cv_res, time_res)
ghs_results <- list(precision = precision, 
                    n = n_res,
                    osb = bh_res)

save(ghs_results, file="ghs_results_stratvns.Rdata")

###########################TB5 5 STRATA#########################################

# data prep
strata_num <- 5
y <- as.vector(tb5$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
min_tb5 <- as.numeric(min_cv[2])

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  res <- STRATVNS(tb5$y, L=strata_num, cvt=0.00336723, imax=5, maxstart=2, tmax=5, kmax = 3, nIterWithNoImpMax=2, nsols = 10, parallelize=TRUE, cputime=1000)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  
  nh <- as.vector(res$nh)
  n <- sum(nh)
  OSB <- as.vector(res$bk)
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
  time_res <- rbind(time_res, time)
  
}


precision <- cbind(var_res, cv_res, time_res)
tb5_results <- list(precision = precision, 
                    n = n_res,
                    osb = bh_res)

save(tb5_results, file="tb5_results_stratvns.Rdata")

###########################GHS 2023 5 STRATA########################################

# data prep
strata_num <- 5
y <- as.vector(ghs23$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
min_ghs23 <- as.numeric(min_cv[4])

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  res <- STRATVNS(ghs23$y, L=strata_num, cvt=0.004899189, imax=5, tmax=5, kmax = 3, maxstart=2, nIterWithNoImpMax=2, nsols = 10, parallelize=TRUE, cputime=1000)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  
  nh <- as.vector(res$nh)
  n <- sum(nh)
  OSB <- as.vector(res$bk)
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
  time_res <- rbind(time_res, time)
  
}


precision <- cbind(var_res, cv_res, time_res)
ghs23_results <- list(precision = precision, 
                      n = n_res,
                      osb = bh_res)

save(ghs23_results, file="ghs23_results_stratvns.Rdata")
