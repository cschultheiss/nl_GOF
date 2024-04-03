fisher.split <- function(sels){
  # function to calculate p-values from Fisher proportion test
  # Input
  # sels (1/0, matrix): 1 if not chosen by FOCI on a given split
  # Output
  # pv (numeric, vector): p-values from proportion test
  B <- nrow(sels) # number of splits
  pv <- rep(1, ncol(sels))
  mn <- mean(sels) # non-selection rate
  cs <- colSums(sels) # sum of non-selection per variable
  rej <- which(cs <= (mn * B + 1e-5)) # reject over-average selected variables
  comp <- max(cs[cs <= (mn * B + 1e-5)]) # value to compare to
  if(any(cs > (mn * B + 1e-5))){
    # proportion test
    pv[-rej] <- sapply(cs[-rej],
                       function(x) fisher.test(cbind(c(B - x, x), c(B - comp, comp)),
                                               alternative = "less")$p.value)
  }
  pv
}

wilc.split <- function(sels){
  # function to calculate p-values from Fisher proportion test
  # Input
  # sels (1/0, matrix): 1 if not chosen by FOCI on a given split
  # Output
  # pv (numeric, vector): p-values from proportion test
  pv <- rep(1, ncol(sels))
  mn <- mean(sels) # non-selection rate
  cm <- colMeans(sels) # average of non-selection per variable
  rej <- which(cm <= mn + 1e-6) # reject over-average selected variables
  max.rej <- which.max(cm[rej]) # variable to compare to
  if(any(cm > (mn + 1e-6))){
    # proportion test
    pv[-rej] <- apply(matrix(sels[,-rej], nrow(sels)), 2,
                      function(x) wilcox.test(x, sels[,rej[max.rej]], alternative = "greater")$p.value)
  }
  pv
}