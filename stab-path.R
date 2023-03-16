require(latex2exp)
require(modeest)
folder <- "results/16-Mar-2023 13.42"
flz <- list.files(folder)
n.lim <- 200
s <- 0
lims.list <- list()
all.s.list <- list()
ns <- numeric(length(flz))
for (file in flz){
  s <- s + 1
  load(paste(folder, "/", file, sep = ""))
  p <- dim(simulation$steps.in)[2]
  
  steps.in <- simulation$steps.in
  steps.out <- simulation$steps.out
  
  sel.in <- simulation$sel.in
  sel.out <- simulation$sel.out
  
  lims <- c(0, unique(sort(c(simulation$steps.in[!is.na(simulation$steps.in)],
                             simulation$steps.out[!is.na(simulation$steps.out)]))))

  if (length(lims > n.lim)){
    lims <- seq(0, max(lims), length.out = n.lim)
  }

  
  steps.in.sorted <- array(NA, dim(steps.in))
  steps.out.sorted <- array(NA, dim(steps.in))
  
  for (i in 1:dim(steps.in)[1]){
    for (j in 1:dim(steps.in)[3]){
      steps.in.sorted[i,,j] <- steps.in[i,sel.in[i,,j],j]
      steps.out.sorted[i,,j] <- steps.out[i,sel.out[i,,j],j]
    }
  }
  reset <- function(steps, lim){
    wi <- which(steps < lim)
    if(length(wi) > 0){
      steps[min(wi):length(steps)] <- NA
    }
    return(steps)
  }
  
  all.s.in <- all.s.out <- array(NA, dim = c(dim(steps.in), length(lims)))
  all.s.in[] <- sapply(lims, function(lim) aperm(apply(steps.in.sorted, c(1,3), reset, lim = lim), c(2,1,3)))
  all.s.out[] <- sapply(lims, function(lim) aperm(apply(steps.out.sorted, c(1,3), reset, lim = lim), c(2,1,3)))
  all.s.var.in <- all.s.var.out <- array(NA, dim(all.s.in))
  for (i in 1:dim(steps.in)[1]){
    for (j in 1:dim(steps.in)[3]){
      sel.in.ij <- sel.in[i,,j]
      sel.in.ij[is.na(sel.in.ij)] <- setdiff(1:p, sel.in.ij)
      all.s.var.in[i,,j,] <- all.s.in[i, order(sel.in.ij),j,]
      
      sel.out.ij <- sel.out[i,,j]
      sel.out.ij[is.na(sel.out.ij)] <- setdiff(1:p, sel.out.ij)
      all.s.var.out[i,,j,] <- all.s.out[i, order(sel.out.ij),j,]
    }
  }
  lims.list[[s]] <- lims
  all.s.list[["in"]][[s]] <- all.s.var.in
  all.s.list[["out"]][[s]] <- all.s.var.out
  ns[s] <- simulation$n
  
  rm(all.s.in, all.s.out, all.s.var.in, all.s.var.out, sel.in.ij, sel.out.ij, steps.in.sorted,
     steps.out.sorted, steps.in, steps.out, sel.in, sel.out, lims)
}


par(mfrow = c(2,2))
for (s in 1:length(flz)){
  lims <- lims.list[[s]]
  add = FALSE
  for (split in c("out", "in")){
    all.s.var <- all.s.list[[split]][[s]]
    frac.run <- apply(!is.na(all.s.var), 2:4, mean)
    frac.loc <- apply(frac.run, c(1,3), mean)
    matplot(lims, t(frac.loc), type = "l", main = ns[s], lty = 1 + 1 * add, ylab = "selection fraction",
            xlab = TeX("$\\delta$ limit"), ylim = c(0, 1), add = add)
    # frac.bounds <- apply(frac.run, c(1,3), quantile, probs = c(0.025, 0.975))
    # for(j in 1:p){
    #   matplot(lims, t(frac.bounds[,j,]), type = "l", main = simulation$n, lty = 2, col = j, add = TRUE)
    # }
    add = TRUE
  }

  # abline(h = 0.5, lty = 2, col = "gray")
  # abline(v = 0.052, lty = 2, col = "gray")
}

# par(mfrow = c(2,2))
# for (file in flz){
#   load(paste(folder, "/", file, sep = ""))
#   steps.in <- simulation$steps.in
#   hist(steps.in[, 1,])
# }

par(mfrow = c(2,2))
for (file in flz){
  load(paste(folder, "/", file, sep = ""))
  steps.in <- simulation$steps.in
  sel.in <- simulation$sel.in
  steps.in.s <- steps.in[, 1,][sel.in[,1,] == 1]
  hist(steps.in.s, main = mlv(steps.in.s, method = "meanshift", na.rm = TRUE)[1], freq = FALSE)
}