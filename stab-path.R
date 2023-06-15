require(latex2exp)
require(modeest)
require(scales)
source("split.R")

folder <- "results/gam-sd0"
savefolder <- "Figures/two-branchf"
flz <- list.files(folder)
nf <- length(flz)
analysis <- paste(folder, "/analysis.RData", sep = "")
test.in <- FALSE
get.pval <- TRUE
if (file.exists(analysis)){
  nf <- nf - 1
  load(analysis)
} else {
  ns <- numeric(nf)
  n.lim <- 200
  s <- 0
  lims.list <- list()
  all.s.list <- list()
  pval.list <- list()
  for (file in flz){
    s <- s + 1
    load(paste(folder, "/", file, sep = ""))
    p <- dim(simulation$steps.out)[2]
    
    if(get.pval) pval.list[[s]] <- simulation$pval.corr
    
    steps.in <- simulation$steps.in
    steps.out <- simulation$steps.out
    
    sel.in <- simulation$sel.in
    sel.out <- simulation$sel.out
    
    lims <- c(0, unique(sort(c(simulation$steps.in[!is.na(simulation$steps.in)],
                               simulation$steps.out[!is.na(simulation$steps.out)]))))
    
    if (length(lims > n.lim)){
      lims <- seq(0, max(lims), length.out = n.lim)
    }
    
    
    steps.in.sorted <- array(NA, dim(steps.out))
    steps.out.sorted <- array(NA, dim(steps.out))
    
    for (i in 1:dim(steps.out)[1]){
      for (j in 1:dim(steps.out)[3]){
        if(test.in) steps.in.sorted[i,,j] <- steps.in[i,sel.in[i,,j],j]
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
    
    all.s.in <- all.s.out <- array(NA, dim = c(dim(steps.out), length(lims)))
    if(test.in) all.s.in[] <- sapply(lims, function(lim) aperm(apply(steps.in.sorted, c(1,3), reset, lim = lim), c(2,1,3)))
    all.s.out[] <- sapply(lims, function(lim) aperm(apply(steps.out.sorted, c(1,3), reset, lim = lim), c(2,1,3)))
    all.s.var.in <- all.s.var.out <- array(NA, dim(all.s.in))
    for (i in 1:dim(steps.out)[1]){
      for (j in 1:dim(steps.out)[3]){
        if(test.in) {
          sel.in.ij <- sel.in[i,,j]
          sel.in.ij[is.na(sel.in.ij)] <- setdiff(1:p, sel.in.ij)
          all.s.var.in[i,,j,] <- all.s.in[i, order(sel.in.ij),j,]
        }
        
        sel.out.ij <- sel.out[i,,j]
        sel.out.ij[is.na(sel.out.ij)] <- setdiff(1:p, sel.out.ij)
        all.s.var.out[i,,j,] <- all.s.out[i, order(sel.out.ij),j,]
      }
    }
    lims.list[[s]] <- lims
    all.s.list[["in"]][[s]] <- all.s.var.in
    all.s.list[["out"]][[s]] <- all.s.var.out
    ns[s] <- simulation$n
    
    rm(all.s.in, all.s.out, all.s.var.in, all.s.var.out, sel.out.ij, steps.in.sorted,
       steps.out.sorted, steps.in, steps.out, sel.in, sel.out, lims)
    if(test.in) rm(sel.in.ij)
  }
  if(get.pval){
    save(all.s.list, lims.list, ns, pval.list, file = analysis)
  } else {
    save(all.s.list, lims.list, ns, file = analysis)
  }
  
}


reverse <- FALSE
sim.sel <- FALSE
add.split <- TRUE
B <- dim(all.s.list$out[[1]])[1]
p <- dim(all.s.list$out[[1]])[2]
stab <- c(1,5)
unstab <- (1:p)[-stab]

par(mfrow = c(2,2))
for (s in 1:nf){
  lims <- lims.list[[s]]
  add = FALSE
  splits <- "out"
  if (test.in) splits <- c("out", "in")
  for (split in splits){
    all.s.var <- all.s.list[[split]][[s]]
    if (sim.sel) {
      if (reverse){
        all.s.var <- pmax(all.s.var[1:(B/2),,,], all.s.var[(B/2 + 1):B,,,], na.rm = TRUE)
      } else {
        all.s.var <- all.s.var[1:(B/2),,,] * all.s.var[(B/2 + 1):B,,,]
      }
      
    }
    frac.run <- apply(!is.na(all.s.var), 2:4, mean)
    frac.loc <- apply(frac.run, c(1,3), mean)
    frac.stab <- apply(matrix(frac.loc[stab,], nrow = length(stab)), 2, mean)
    frac.unstab <- apply(matrix(frac.loc[unstab,], nrow = length(unstab)), 2, mean)
    # matplot(lims, t(frac.loc), type = "l", main = ns[s], lty = 1 + 1 * add, ylab = "selection fraction",
    #         xlab = TeX("$\\delta$ limit"), ylim = c(0, 1), add = add)
    matplot(lims, cbind(frac.stab, frac.unstab), type = "l", main = ns[s], lty = 1 + 1 * add, ylab = "selection fraction",
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

# ecdf
par(mfrow = c(2,2))
for (s in 1:nf){
  all.s.0 <- all.s.list$out[[s]][,,,1]
  if (sim.sel) {
    if(reverse) {
      all.s.0 <- pmax(all.s.0[1:(B/2),,], all.s.0[(B/2 + 1):B,,], na.rm = TRUE)
    } else {
      all.s.0 <- all.s.0[1:(B/2),,] * all.s.0[(B/2 + 1):B,,]
    }
    
  }
  
  all.frac <- apply(!is.na(all.s.0), 2:3, mean)  
  plot.ecdf(all.frac[stab,], xlim = c(0, 1), main = ns[s])
  if (length(unstab) > 0)
    plot.ecdf(all.frac[unstab,], col = 2, add= TRUE)
    # abline(v =0.75)
    # abline(h = 0.5)
}

#ROC
par(mfrow = c(2,2))
for (s in 1:nf){
  all.s.0 <- all.s.list$out[[s]][,,,1]
  if (get.pval){
    pval.lim <- 0.05
    pval <- pval.list[[s]]
    glob <- which(pval > pval.lim)
  }
    
  
  Bb <- B
  pi0 <- 0.5
  if (sim.sel) {
    if(reverse) {
      all.s.0 <- pmax(all.s.0[1:(B/2),,] , all.s.0[(B/2 + 1):B,,], na.rm = TRUE)
      pi0 <- 0.75
    } else {
      all.s.0 <- all.s.0[1:(B/2),,] * all.s.0[(B/2 + 1):B,,]
      pi0 <- 0.25
    }
    Bb <- Bb / 2
  }

  pis <- sort(unique(c((0:Bb)/Bb, pi0)))

  avg.frac <- apply(is.na(all.s.0), 3, mean)
  all.frac <- apply(is.na(all.s.0), 2:3, mean)
  if (get.pval) all.frac[, glob] <- 1
  r1 <- sapply(pis, function(pi) mean(all.frac[unstab,] >= pi))
  r2 <- sapply(pis, function(pi) mean(all.frac[stab,] >= pi))
  if(length(unstab) < 1)
    cat("Found ", 100 * r2[pis == pi0], "% of the stable variables with the 50 % requirement \n", sep ="")
  # bds <- sapply(round(avg.frac, 2), function(avg) bound(avg, B, p, 0.1))
  which.sel <- t(all.frac) > (avg.frac)
  if(length(unstab) < 1)
    cat("Found ", 100 * mean(which.sel), "% of the stable variables by comparing to the mean \n", sep ="")
  fdr <- sapply(pis, function(pi) mean(apply(matrix(all.frac[unstab,], nrow = length(unstab)) >= pi, 2, sum) /
                                         pmax(apply(all.frac >= pi, 2, sum), 1)))
  # r1 <- fdr
  # which.sel.tau <- t(all.frac) > bds
  plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", main = ns[s])
  points(r1[pis == pi0], r2[pis == pi0], pch = 4)
  points(mean(is.na(all.s.list$out[[s]][1,unstab,,1])), mean(is.na(all.s.list$out[[s]][1,stab,,1])), pch = 2, col = 4)
  points(mean(which.sel[,unstab]), mean(which.sel[,stab]), col = 3, pch = 3)
  # points(mean(which.sel.tau[,unstab]), mean(which.sel.tau[,stab]), col = 4, pch = 4)
  abline(0,1, col = "gray", lty= 2)
  
  if(add.split){
    pv <- apply(1 * is.na(all.s.0), 3, fisher.split)
    pi0 <- 1e-2
    pis <- sort(unique(c(pv, pi0)))
    if(get.pval) pv[, glob] <- 0
    if(length(unstab) < 1)
      cat("Found ", 100 * mean(pv <= pi0), "% of the stable variables with the splitting test \n", sep ="")
    fdr <- sapply(pis, function(pi) mean(apply(matrix(pv[unstab,], nrow = length(unstab)) <= pi, 2, sum) / pmax(apply(pv <= pi, 2, sum), 1)))
    r1 <- sapply(pis, function(pi) mean(pv[unstab,] <= pi))
    r2 <- sapply(pis, function(pi) mean(pv[stab,] <= pi))
    # r1 <- fdr
    lines(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", main = ns[s], col = 2)
    points(r1[pis == pi0], r2[pis == pi0], pch = 4, col = 2)
  }
}

#ROC alt
# png(paste(savefolder, "/ROC+ECDF+box.png", sep = ""), width = 3600,
# height = 1200, res = 300)
par(mfrow = c(1, 3))
exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
for (s in 1:nf){
  all.s.0 <- all.s.list$out[[s]][,,,1]
  Bb <- B
  pi0 <- 1e-2
  
  if (get.pval){
    pval.lim <- 0.05
    pval <- pval.list[[s]]
    glob <- which(pval > pval.lim)
    load(paste(folder, "/", flz[s+1], sep = ""))
    glob.s <- t(simulation$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
  }
  if (sim.sel) {
    if(reverse) {
      all.s.0 <- pmax(all.s.0[1:(B/2),,] , all.s.0[(B/2 + 1):B,,], na.rm = TRUE)
      pi0 <- 0.75
    } else {
      all.s.0 <- all.s.0[1:(B/2),,] * all.s.0[(B/2 + 1):B,,]
      pi0 <- 0.25
    }
    Bb <- Bb / 2
  }

  pv <- apply(1 * is.na(all.s.0), 3, fisher.split)
  if(get.pval) pv[, glob] <- 0
  pis <- sort(unique(c(pv)))
  
  r1 <- sapply(pis, function(pi) mean(pv[unstab,] <= pi))
  r2 <- sapply(pis, function(pi) mean(pv[stab,] <= pi))
  
  if (s == 1) {
    plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[s], lty = ltys[s],
         xlab = "False positive rate", ylab = "True positive rate", cex.lab = exp.text, lwd = exp.lines)
  } else {
    points(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[s], lty = ltys[s], lwd = exp.lines)
  }
  
  points(mean(pv[unstab,] <= pi0), mean(pv[stab,] <= pi0), pch = 4, col = cols[s], cex = exp.points)
  # points(mean(is.na(all.s.list$out[[s]][1,unstab,,1])), mean(is.na(all.s.list$out[[s]][1,stab,,1])), pch = 2, col = 2)
  t1.mat <- is.na(all.s.list$out[[s]][,unstab,,1])
  p.mat <- is.na(all.s.list$out[[s]][,stab,,1])
  if (get.pval){
    if(length(dim(t1.mat)) > 2)
      t1.mat <- aperm(t1.mat, c(1, 3, 2))
    if(length(dim(p.mat)) > 2)
      p.mat <- aperm(p.mat, c(1, 3, 2))
    t1.mat <- pmax(t1.mat, glob.s)
    p.mat <- pmax(p.mat, glob.s)
  }


  points(mean(t1.mat), mean(p.mat), pch = pchs[s], col = cols[s], cex = exp.points)
  # points(mean(which.sel.tau[,unstab]), mean(which.sel.tau[,stab]), col = 4, pch = 4)
  abline(0,1, col = "gray", lty= 2)
  
  labels <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns), "$')", sep = "", collapse = ","), ")")))
  legend('bottomright', legend = labels, col = cols, lty = ltys, pch = pchs, cex = exp.text, pt.cex = 1, lwd = exp.lines)
}

# par(mfrow = c(1,1))
cols <- 1:4
ltys <- 1:4
s <- 0
ns.p <- numeric(0)
for (file in flz){
  if(grepl("results", file)){
    s <- s + 1
    load(paste(folder, "/", file, sep = ""))
    if(max(simulation$pval) < 1e-3) next()
    ns.p <- c(ns.p, simulation$n)
    np <- prod(dim(simulation$pval))
    npc <- length(simulation$pval.corr)
    if (simulation$n == min(ns)){
      plot(c(sort(simulation$pval), 1), (1:(np + 1))/(np + 1), type = "l", xlim = c(0, 1),
           col = cols[s], lty = ltys[s], xlab = "p", ylab = "Fn(p)", cex.lab = exp.text, lwd = exp.lines)
    } else {
      lines(c(sort(simulation$pval), 1), (1:(np + 1))/(np + 1), col = cols[s], lty = ltys[s], lwd = exp.lines)
    }
    
    lines(c(sort(simulation$pval.corr), 1), (1:(npc + 1))/(npc + 1), col = cols[s], lty = ltys[s], lwd = exp.lines)
    points(c(sort(simulation$pval.corr), 1), (1:(npc + 1))/(npc + 1), col = alpha(cols[s], 0.2))
    # plot.ecdf(simulation$pval, xlim = c(0,1), main = paste("10^", log10(simulation$n), sep = ""))
    # plot.ecdf(simulation$pval.corr, col = 2, add = TRUE)
  }
}
labels.sub <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns.p), "$')", sep = "", collapse = ","), ")")))
legend('bottomright', legend = labels.sub, col = cols, lty = ltys, cex = exp.text, pt.cex = 1, lwd = exp.lines)
# dev.off()

mspe <- 1
s <- 0
mses <- list()
rcor <- list()
for (file in flz){
  if(grepl("results", file)){
    s <- s + 1
    load(paste(folder, "/", file, sep = ""))
    mses[[s]] <- simulation$all[,1]
    rcor[[s]] <- simulation$all[,2]
  }
}
mse <- matrix(unlist(mses), ncol = nf)
boxplot(mse/mspe, log = "y", names = labels, ylab = "Relative approximation error", cex.lab = exp.text, cex.axis =exp.text, yaxt ="n")
axis(side = 2)
# dev.off()


# only with two predictors
mn <- Vectorize(function(i, pi = pi0){
  sr <- TRUE
  j <- i - 1
  while (j >= 0 && sr) {
    if (fisher.split(cbind(c(rep(1, i), rep(0, B-i)), c(rep(1, j), rep(0, B-j))))[1] < pi){
      sr <- FALSE
    } else {
      j <- j - 1
    }
  }
 j
})

par(mfrow = c(2,2))
for (s in 1:nf){
  all.s.0 <- all.s.list$out[[s]][,,,1]
  ina <- apply(is.na(all.s.0), 2:3, sum)
  plot(ina[unstab,], ina[stab,], xlim = c(0, B), ylim = c(0, B),
       col = 1 + 2 *(ina[unstab,] <= mn(ina[stab,])) + 1 *(ina[stab,] <= mn(ina[unstab,])),
       xlab = "unstable", ylab = "stable")
  lines(0:B, mn(0:B))
  lines(mn(0:B), 0:B)
}

# par(mfrow = c(2,2))
# for (file in flz){
#   load(paste(folder, "/", file, sep = ""))
#   steps.in <- simulation$steps.in
#   sel.in <- simulation$sel.in
#   steps.in.s <- steps.in[, unstab,][sel.in[,1,] %in% unstab]
#   hist(steps.in.s, main = mlv(steps.in.s, method = "meanshift", na.rm = TRUE)[1], freq = FALSE)
# }