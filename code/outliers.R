detectOutliers <- function(x, t, nPC) {
  pca <- prcomp(t(na.omit(x)))
  outliers <- foreach ( pc = 1:nPC) %do% {
      abs(pca$x[,pc] - mean(pca$x[, pc])) > t * sd(pca$x[, pc])
    } %>%
    do.call("rbind", .) %>% 
    apply(2, any)     
  return(list(pca=pca, outliers=outliers))
}

plotPCAPairs <- function(pca, nPC = 3, color = NULL, shape = NULL,
  main = "Principal component analysis") 
{
  if (length(color) == 0) color <- "black"
  if (length(shape) == 0) shape <- 1
  pairs(pca$x[, 1:nPC], 
    col=as.factor(color),
    pch = as.factor(shape) %>% as.numeric, 
    main=main)
}

removeOutliers <- function(data, outliers) {
  shrinkList(data, !outliers, nrow(data$sample))
}



