## Analysis of GA4Stratification package for computational methods in SA context ##
# Date: 08 June 2025
# Author: Georgi Borros

# Set working directory to path with 'reproducible-code' folder
setwd("C:/Users/01459189/OneDrive/phd/r/sa_comp_methods/reproducible-code")

library(sampling)
library(data.table)
library(stratbr)
library(haven)
library(ggplot2)
options(scipen=999)

# Source author-written functions
source("dfstrata.R")
source("elbow_plot.R")

################################################################################
## Part 1: Load Datasets
################################################################################

##### Dataset 1: ECD Census 2021 ######
load("ecd_census.RData")
df <- df[order(df$fees_amount),]
ecd <- as.data.frame(df[, c("fees_amount")])
names(ecd) <- "y"
ecd <- as.data.frame(ecd[ecd$y!=0,])
names(ecd) <- "y"
ecd_n <- 3711 ## 10% sample size, rounded up

#### Dataset 2: GHS 2024 ####

load("ghs_2024.RData")
df <- df[order(df$totmhinc),]
ghs <- as.data.frame(df[, c("totmhinc")])
names(ghs) <- "y"
ghs_n <- 2051


#### Dataset 3: GHS 2023 ####

load("ghs_2023.RData")
df <- df[order(df$totmhinc),]
ghs23 <- as.data.frame(df[, c("totmhinc")])
names(ghs23) <- "y"
ghs23_n <- 2055
set.seed(1234)
elbow_plot(ghs23, 8)


#### Dataset 4: TB5 ####
load("tb5.RData")
df <- df[df$age<=59,]
df <- df[order(df$total_elom),]
tb5<- as.data.frame(df[, c("total_elom")])
names(tb5) <- "y"
tb5_n <- 523

## Change wd to GA4Stratification folder (saved within reproducible-code folder)
setwd("C:/Users/01459189/OneDrive/phd/r/sa_comp_methods/reproducible-code/GA4Stratification")
source("GA4Stratification_updated.R")

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

# algorithm
for (i in 1:10) {
  set.seed(i)
start.time <- Sys.time()
res <- GA4Stratification(ecd, ## data
                         5, ## strata
                         ecd_n, ## sample size
                         50, ## pop size
                         5000, ## iterations
                         0.2, ## mutation rate
                         4) ## allocation method
end.time <- Sys.time()
time <- difftime(end.time, start.time, units="mins")

nh <- as.vector(res[,2])
## Boundaries
OSB_1 <- as.vector(res[,1])
OSB <- numeric()
OSB[1] <- OSB_1[1]
for (l in 2:strata_num) {
  
  OSB[l] <- OSB[l-1] + OSB_1[l]
  
  print(OSB)
}
bh <- numeric()

for (l in 1:(strata_num-1)) {
  bh[l] <- y[OSB[l]]
}
# Put results into the dataframe
df <- dfstrata(df=ecd, X=y, h= strata_num, OSB=bh, nh=nh)

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

save(ecd_results, file="C:/Users/SALDRU-SDH/Documents/GitHub/sa_comp_methods/ecd_results_GA4.Rdata")

###########################GHS 5 STRATA########################################

# data prep
strata_num <- 5
y <- as.vector(ghs$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()

# algorithm
# algorithm
for (i in 1:10) {
  set.seed(i)
  start.time <- Sys.time()
  res <- GA4Stratification(ghs, ## data
                           5, ## strata
                           ghs_n, ## sample size
                           50, ## pop size
                           5000, ## iterations
                           0.2, ## mutation rate
                           4) ## allocation method
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")

  nh <- as.vector(res[,2])
  ## Boundaries
  OSB_1 <- as.vector(res[,1])
  OSB <- numeric()
  OSB[1] <- OSB_1[1]
  for (l in 2:strata_num) {
    
    OSB[l] <- OSB[l-1] + OSB_1[l]
    
    print(OSB)
  }
  bh <- numeric()
  
  for (l in 1:(strata_num-1)) {
    bh[l] <- y[OSB[l]]
  }
# Put results into the dataframe
df <- dfstrata(df=ghs, X=y, h= strata_num, OSB=bh, nh=nh)

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

save(ghs_results, file="C:/Users/SALDRU-SDH/Documents/GitHub/sa_comp_methods/ghs_results_GA4.Rdata")

#############################TB5 5 STRATA######################################

# data prep
strata_num <- 5
y <- as.vector(tb5$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()

# algorithm
for (i in 1:10) {
  set.seed(i)
start.time <- Sys.time()
res <- GA4Stratification(tb5, ## data
                         5, ## strata
                         tb5_n, ## sample size
                         50, ## pop size
                         5000, ## iterations
                         0.2, ## mutation rate
                         4) ## allocation method
end.time <- Sys.time()
time <- difftime(end.time, start.time, units="mins")

nh <- as.vector(res[,2])
## Boundaries
OSB_1 <- as.vector(res[,1])
OSB <- numeric()
OSB[1] <- OSB_1[1]
for (l in 2:strata_num) {
  
  OSB[l] <- OSB[l-1] + OSB_1[l]
  
  print(OSB)
}
bh <- numeric()

for (l in 1:(strata_num-1)) {
  bh[l] <- y[OSB[l]]
}
# Put results into the dataframe
df <- dfstrata(df=tb5, X=y, h= strata_num, OSB=bh, nh=nh)

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

save(tb5_results, file="C:/Users/SALDRU-SDH/Documents/GitHub/sa_comp_methods/tb5_results_GA4.Rdata")
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

# algorithm
# algorithm
for (i in 1:10) {
  set.seed(i)
  start.time <- Sys.time()
  res <- GA4Stratification(ghs23, ## data
                           5, ## strata
                           ghs23_n, ## sample size
                           50, ## pop size
                           5000, ## iterations
                           0.2, ## mutation rate
                           4) ## allocation method
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  
  nh <- as.vector(res[,2])
  ## Boundaries
  OSB_1 <- as.vector(res[,1])
  OSB <- numeric()
  OSB[1] <- OSB_1[1]
  for (l in 2:strata_num) {
    
    OSB[l] <- OSB[l-1] + OSB_1[l]
    
    print(OSB)
  }
  bh <- numeric()
  
  for (l in 1:(strata_num-1)) {
    bh[l] <- y[OSB[l]]
  }
  # Put results into the dataframe
  df <- dfstrata(df=ghs23, X=y, h= strata_num, OSB=bh, nh=nh)
  
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

save(ghs23_results, file="C:/Users/SALDRU-SDH/Documents/GitHub/sa_comp_methods/ghs23_results_GA4.Rdata")



















