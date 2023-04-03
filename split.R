wilc.split <- function(sels){
  pv <- rep(1, ncol(sels))
  mn <- mean(sels)
  cm <- colMeans(sels)
  rej <- which(cm <= mn)
  max.rej <- which.max(cm[rej])
  pv[-rej] <- apply(matrix(sels[,-rej], nrow(sels)), 2,
        function(x) wilcox.test(x, sels[,rej[max.rej]], alternative = "greater")$p.value)
  pv
}