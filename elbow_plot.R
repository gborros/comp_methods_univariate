elbow_plot <- function(df, k) {
  result <- numeric()
  
  for (i in 1:k){
    
    km <- kmeans(df, i, nstart=50)
    twss <- km$tot.withinss
    result <- rbind(result, twss)
    
  }
  
  clusters <- c(1:k)

  p <- as.data.frame(result,clusters)
  ggplot(data=p, aes(x=clusters, y=V1, group=1)) +
    geom_line(color="gray56", size=1)+ scale_x_continuous(breaks =c(1:k)) +
    geom_point() + theme_classic(base_family="serif", base_size=8) + labs(title="",x="Clusters (k)", y = "") + 
    geom_vline(xintercept = 5, linetype="dashed", 
                    color = "palevioletred", size=1)
}