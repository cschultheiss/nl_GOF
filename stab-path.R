require(latex2exp)
folder <- "results/17-Feb-2023 17.38"
flz <- list.files(folder)
s <- 0
lims.list <- list()
all.s.list <- list()
ns <- numeric(length(flz))
for (file in flz){
  s <- s + 1
  load(paste(folder, "/", file, sep = ""))
  p <- dim(simulation$steps)[2]
  
  steps <- simulation$steps
  
  sel <- simulation$sel
  
  lims <- c(0, unique(sort(simulation$steps[!is.na(simulation$steps)])))

  if (length(lims > 1000)){
    lims <- seq(0, max(lims), length.out = 1000)
  }

  
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
  
  all.s <- array(NA, dim = c(dim(steps), length(lims)))
  all.s[] <- sapply(lims, function(lim) aperm(apply(steps.sorted, c(1,3), reset, lim = lim), c(2,1,3)))
  all.s.var <- array(NA, dim(all.s))
  for (i in 1:dim(steps)[1]){
    for (j in 1:dim(steps)[3]){
      sel.ij <- sel[i,,j]
      sel.ij[is.na(sel.ij)] <- setdiff(1:p, sel.ij)
      all.s.var[i,,j,] <- all.s[i, order(sel.ij),j,]
    }
  }
  lims.list[[s]] <- lims
  all.s.list[[s]] <- all.s.var
  ns[s] <- simulation$n
}

par(mfrow = c(2,2))
for (s in 1:length(flz)){
  lims <- lims.list[[s]]
  all.s.var <- all.s.list[[s]]
  frac.run <- apply(!is.na(all.s.var), 2:4, mean)
  frac.loc <- apply(frac.run, c(1,3), mean)
  matplot(lims, t(frac.loc), type = "l", main = ns[s], lty = 1, ylab = "Selection fraction",
          xlab = TeX("$\\delta$ limit"), ylim = c(0, 1))
  frac.bounds <- apply(frac.run, c(1,3), quantile, probs = c(0.025, 0.975))
  for(j in 1:p){
    matplot(lims, t(frac.bounds[,j,]), type = "l", main = simulation$n, lty = 2, col = j, add = TRUE)
  }
}

par(mfrow = c(2,2))
for (file in flz){
  load(paste(folder, "/", file, sep = ""))
  steps <- simulation$steps
  hist(steps[, 1,])
}