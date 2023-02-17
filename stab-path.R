steps <- simulation$steps

sel <- simulation$sel

lims <- c(0, unique(sort(simulation$steps[!is.na(simulation$steps)])))

steps.sorted <- array(NA, dim(steps))

for (i in 1:dim(steps)[1]){
  for (j in 1:dim(steps)[3]){
    steps.sorted[i,,j] <- steps[i,sel[i,,j],j]
  }
}

reset <- function(steps, lim){
  wi <- which(steps < lim)
  if(length(wi) > 0){
    steps[min(wi):length(steps)] <- NA
  }
  return(steps)
}
aperm(apply(steps.sorted, c(1,3), reset, lim = 0), c(2,1,3))

all.s <- array(NA, dim = c(dim(steps), length(lims)))
all.s[] <- sapply(lims, function(lim) aperm(apply(steps.sorted, c(1,3), reset, lim = lim), c(2,1,3)))
all.s.var <- array(NA, dim(all.s))
for (i in 1:dim(steps)[1]){
  for (j in 1:dim(steps)[3]){
    all.s.var[i,,j,] <- all.s[i, order(sel[i,,j]),j,]
  }
}
frac.run <- apply(!is.na(all.s.var), 2:4, mean)
frac.fun <- apply(frac.run, c(1,3), mean)
matplot(lims, t(frac.fun), type = "l")
