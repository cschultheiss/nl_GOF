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

fisher.split <- function(sels){
  B <- nrow(sels)
  pv <- rep(1, ncol(sels))
  mn <- mean(sels)
  cs <- colSums(sels)
  rej <- which(cs <= mn * B)
  comp <- max(cs[cs <= mn * B])
  if(any(cs > mn * B))
    pv[-rej] <- sapply(cs[-rej],
                    function(x) fisher.test(cbind(c(B - x, x), c(B - comp, comp)),
                                                  alternative = "less")$p.value)
  pv
}