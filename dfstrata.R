###----------------### function to put strata in dataframe ###-----------------

dfstrata <- function(df, ## dataframe
                     X, ## stratification variable
                     h, ## number of strata
                     OSB, ## OSB vector
                     nh) ## nh vector
{
  
  strata <- c("strata1")
  df$strata <- strata
  df$X <- X
  
  for (i in 1:(h-1)) {
    df$strata[df$X > OSB[i]] <- paste0("strata",i+1)
  }
  
  n <- c(1)
  df$n <- n
  
  for (i in 1:h) {
    df$n[df$strata == paste0("strata",i)] <- nh[i]
  }
  
  df <- df[order(df$strata,df$strata),]
  return(df)
}