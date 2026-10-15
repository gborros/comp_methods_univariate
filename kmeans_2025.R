## Analysis of k-means method for stratification in SA context ##
# Date: 18 July 2025
# Author: Georgi Borros

# Set working directory to path with 'reproducible-code' folder
setwd("C:/Users/01459189/OneDrive/phd/r/sa_comp_methods/reproducible-code")

library(sampling)
library(data.table)
library(stratification)
library(factoextra)
library(NbClust)
library(haven)
library(dplyr)
options(scipen=999)

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
#ecd_n <- 3711 ## 10% sample size, rounded up
ecd_n <- 8400 ## 8400 to acheive cv of 0.0037 (comment this out when looking at sample size of 10% and use the line above)

#### Dataset 2: GHS 2024 ####

load("ghs_2024.RData")
df <- zap_label(zap_labels(df[order(df$totmhinc),]))
ghs <- as.data.frame(df[, c("totmhinc")])
names(ghs) <- "y"
#ghs_n <- 2051
ghs_n <- 10300 ## 10300 to achieve cv of 0.0049


#### Dataset 3: GHS 2023 ####

load("ghs_2023.RData")
df <- zap_label(zap_labels(df[order(df$totmhinc),]))
ghs23 <- as.data.frame(df[, c("totmhinc")])
names(ghs23) <- "y"
#ghs23_n <- 2055
ghs23_n <- 11800 ## 11800 to achieve cv around 0.0049


#### Dataset 4: TB5 ####
load("tb5.RData")
df <- zap_label(zap_labels(df[df$age<=59,]))
df <- df[order(df$total_elom),]
tb5<- as.data.frame(df[, c("total_elom")])
names(tb5) <- "y"
#tb5_n <- 523
tb5_n <- 540 ## 540 needed to achieve cv 0.0034

################################################################################
## Part 2: Generate Results for each dataset
################################################################################

###########################ECD Census 5 STRATA#######################################
strata_num <- 5
y <- as.vector(ecd$y)
N <- length(y)
ssize <- ecd_n
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
strata <- list()

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
km <- kmeans(ecd, 5, nstart=50)
end.time <- Sys.time()
time <- difftime(end.time, start.time, units="mins")
ecd$strata <- km$cluster
df <- ecd
dfg <- tibble()
dfg <- df %>%
  group_by(strata) %>%
  summarise(mean = mean(y), sd = sd(y),
            N_h = n(), .groups = "drop")

N <- dfg$N_h 
W <- dfg$N_h / sum(N) ## population proportion of stratum
L <- nrow(dfg) ## number of strata
S <- dfg %>% select(ends_with("_sd")) %>% as.matrix()

  n <- numeric()
  for (h in 1:L){
    n[h] <- as.integer(ssize*(N[h]/sum(N)))
  }
  
  if (sum(n)>ssize) {
    N=sum(n)
    n[which.min(N-n)] <- n[which.min(N-n)] - (sum(n)-ssize) 
  }
  
  if (sum(n)<ssize) {
    N=sum(n)
    n[which.max(N-n)] <- n[which.max(N-n)] + (ssize-sum(n)) 
  }
  
  for (h in 1:L){
    if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
      n[h] = 2
      n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
    }
  }
  
nh <- n

## Calculate accordingly
count <- aggregate(df$y, list(df$strata), FUN=length)
Nh <- count[,2]
N <- length(y)
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
  select <- df$y[df$strata == l] ## select only obs in the strata
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
time_res <- rbind(time_res, time)

strata[[i]] <- ecd$strata

}

precision <- cbind(var_res, cv_res, time_res)
ecd_results <- list(precision = precision, 
                    n = n_res, strata=strata)
mean(precision[,2])
#save(ecd_results, file="ecd_results_kmeans.Rdata")

###########################GHS 2024 5 STRATA##########################################
strata_num <- 5
y <- as.vector(ghs$y)
N <- length(y)
ssize <- ghs_n
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
strata <- list()

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  km <- kmeans(ghs, 5, nstart=50)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  ghs$strata <- km$cluster
  df <- ghs
  dfg <- tibble()
  dfg <- df %>%
    group_by(strata) %>%
    summarise(mean = mean(y), sd = sd(y),
              N_h = n(), .groups = "drop")
  
  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  S <- dfg %>% select(ends_with("_sd")) %>% as.matrix()
  
  n <- numeric()
  for (h in 1:L){
    n[h] <- as.integer(ssize*(N[h]/sum(N)))
  }
  
  if (sum(n)>ssize) {
    N=sum(n)
    n[which.min(N-n)] <- n[which.min(N-n)] - (sum(n)-ssize) 
  }
  
  if (sum(n)<ssize) {
    N=sum(n)
    n[which.max(N-n)] <- n[which.max(N-n)] + (ssize-sum(n)) 
  }
  
  for (h in 1:L){
    if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
      n[h] = 2
      n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
    }
  }
  
  nh <- n
  
  ## Calculate accordingly
  count <- aggregate(df$y, list(df$strata), FUN=length)
  Nh <- count[,2]
  N <- length(y)
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
    select <- df$y[df$strata == l] ## select only obs in the strata
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
  time_res <- rbind(time_res, time)
  
  strata[[i]] <- ghs$strata
  
}

precision <- cbind(var_res, cv_res, time_res)
ghs_results <- list(precision = precision, 
                    n = n_res,
                    strata=strata)
mean(precision[,2])
#save(ghs_results, file="ghs_results_kmeans.Rdata")

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
ssize <- tb5_n
strata <- list()

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  km <- kmeans(tb5, 5, nstart=50)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  tb5$strata <- km$cluster
  df <- tb5
  dfg <- tibble()
  dfg <- df %>%
    group_by(strata) %>%
    summarise(mean = mean(y), sd = sd(y),
              N_h = n(), .groups = "drop")
  
  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  S <- dfg %>% select(ends_with("_sd")) %>% as.matrix()
  
  n <- numeric()
  for (h in 1:L){
    n[h] <- as.integer(ssize*(N[h]/sum(N)))
  }
  
  if (sum(n)>ssize) {
    N=sum(n)
    n[which.min(N-n)] <- n[which.min(N-n)] - (sum(n)-ssize) 
  }
  
  if (sum(n)<ssize) {
    N=sum(n)
    n[which.max(N-n)] <- n[which.max(N-n)] + (ssize-sum(n)) 
  }
  
  for (h in 1:L){
    if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
      n[h] = 2
      n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
    }
  }
  
  nh <- n
  
  ## Calculate accordingly
  count <- aggregate(df$y, list(df$strata), FUN=length)
  Nh <- count[,2]
  N <- length(y)
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
    select <- df$y[df$strata == l] ## select only obs in the strata
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
  time_res <- rbind(time_res, time)
  
  strata[[i]] <- tb5$strata
  
}

precision <- cbind(var_res, cv_res, time_res)
tb5_results <- list(precision = precision, 
                    n = n_res,
                    strata = strata)
mean(precision[,2])
#save(tb5_results, file="tb5_results_kmeans.Rdata")

###########################GHS 2023 5 STRATA########################################

# data prep
# data prep
strata_num <- 5
y <- as.vector(ghs23$y)
N <- length(y)
var_res <- numeric()
cv_res <- numeric() 
n_res <- numeric()
bh_res <- numeric()
time_res <- numeric()
ssize <- ghs23_n
strata <- list()

for (i in 1:10) {
  # algorithm
  set.seed(i)
  start.time <- Sys.time()
  km <- kmeans(ghs23, 5, nstart=50)
  end.time <- Sys.time()
  time <- difftime(end.time, start.time, units="mins")
  ghs23$strata <- km$cluster
  df <- ghs23
  dfg <- tibble()
  dfg <- df %>%
    group_by(strata) %>%
    summarise(mean = mean(y), sd = sd(y),
              N_h = n(), .groups = "drop")
  
  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  S <- dfg %>% select(ends_with("_sd")) %>% as.matrix()
  
  n <- numeric()
  for (h in 1:L){
    n[h] <- as.integer(ssize*(N[h]/sum(N)))
  }
  
  if (sum(n)>ssize) {
    N=sum(n)
    n[which.min(N-n)] <- n[which.min(N-n)] - (sum(n)-ssize) 
  }
  
  if (sum(n)<ssize) {
    N=sum(n)
    n[which.max(N-n)] <- n[which.max(N-n)] + (ssize-sum(n)) 
  }
  
  for (h in 1:L){
    if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
      n[h] = 2
      n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
    }
  }
  
  nh <- n
  
  ## Calculate accordingly
  count <- aggregate(df$y, list(df$strata), FUN=length)
  Nh <- count[,2]
  N <- length(y)
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
    select <- df$y[df$strata == l] ## select only obs in the strata
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
  time_res <- rbind(time_res, time)
  
  strata[[i]] <- ghs23$strata
  
}


precision <- cbind(var_res, cv_res, time_res)
ghs23_results <- list(precision = precision, 
                      n = n_res,
                      strata=strata)
mean(precision[,2])
#save(ghs23_results, file="ghs23_results_kmeans.Rdata")


