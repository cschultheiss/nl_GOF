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
  if(any(cs > (mn * B + 1e-5))){
    rej <- which(cs <= (mn * B + 1e-5))
    comp <- max(cs[cs <= (mn * B + 1e-5)])
    pv[-rej] <- sapply(cs[-rej],
                       function(x) fisher.test(cbind(c(B - x, x), c(B - comp, comp)),
                                               alternative = "less")$p.value)
  }
  pv
}

bound <- function(theta, B, p, t1){
  mins <- minD(theta, B/2)
  fac <- t1 / p
  ind <- min(which(mins < fac))
  ind / B
}