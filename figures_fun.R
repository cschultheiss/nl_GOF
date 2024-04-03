require(latex2exp)
require(abind)
require(scales)
source("split.R")

figures4 <- function(folder){
  # wrapper function to generate plots from Section 4
  # Input
  # folder (string): location where the result files are stored
  
  # get all results file, i.e., not the setup
  flz <- list.files(folder)
  flz <- flz[grep("results", flz)]
  p.tot <- 5 # number of variables apart from y
  p.mes <- 3 # number of observed covariates
  load(paste(folder, "/", flz[1], sep = ""))
  nsim <- dim(simulation[[4]]$pval)[1] # simulation runs per sample size
  B <- dim(simulation[[4]]$pval)[2] # number of splits used
  
  all.comb <- combn(1:p.tot, p.mes) # all covariate combinations
  
  # find covariates with well-specified effect in each setup
  which.stab <- function(sel){
    # which combination is it
    sel.pos <- which(duplicated(rbind(sel, t(all.comb))))- 1
    if (sel.pos == 1){
      1:3
    } else if(sel.pos %in% c(2, 4, 7)){
      NULL
    } else if(sel.pos == 3){
      c(1, 5)
    } else if(sel.pos == 5){
      c(1, 3, 5)
    } else if(sel.pos %in% c(6, 9, 10)){
      5
    } else if(sel.pos == 8){
      c(3, 5)
    } else {
      error("did not find a match")
    }
  }
  
  # put everything in one vector
  stab.col <- rep(FALSE, prod(dim(all.comb)))
  
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    # find columns with well-specified effect per setting
    stab.col[(s - 1) * p.mes + which(mes %in% which.stab(mes))] <- TRUE
    if(length(which.stab(mes)) == length(mes))
      # treat globally well-specified scenarios separately
      stab.col[(s - 1) * p.mes + 1:p.mes] <- NA
    
  }
  stab <- which(stab.col) # columns with well-specified effect
  unstab <- which(!stab.col) # columns with not well-specified effect
  
  get.pval <- TRUE # apply global test
  
  # plotting parameters
  par(mfrow = c(1, 3))
  exp.text <- 1.5
  exp.points <- 1.5
  exp.lines <- 1.5
  cols <- ltys <-  1:4
  pchs <- c(0:2, 6)
  
  t <- 0
  ns <- numeric(length(flz)) # sample sizes
  for (file in flz){
    # loop through sample sizes
    t <- t + 1
    load(paste(folder, "/", file, sep = ""))
    ns[t] <- simulation$n
    # matrix and array to combine all information
    pv.all <- matrix(NA, ncol = nsim, nrow = 0)
    all.s.all <- array(NA, dim = c(B, 0, nsim))
    for (s in 1:ncol(all.comb)){
      # loop through covariate settings
      mes <- all.comb[, s]
      # FOCI step size for the given setting
      all.s <- simulation[[paste(mes, collapse = "")]]$steps.out
      if (get.pval){
        pval.lim <- 0.05
        # where is global model accepted
        glob <- which(simulation[[paste(mes, collapse = "")]]$pval.corr > pval.lim)
        # if we only consider single splits
        glob.s <- t(simulation[[paste(mes, collapse = "")]]$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
      }
      
      all.s <- is.na(all.s) # value is irrelevant
      pv <- apply(1 * all.s, 3, fisher.split) # p-values from proportion test
      if(get.pval) pv[, glob] <- 0 # set all to 0 if global test is not rejected
      pv.all <- rbind(pv.all, pv) # store to matrix
      
      if (get.pval){
        # consider single splits
        all.s <- aperm(all.s, c(1, 3, 2))
        all.s <- pmax(all.s, glob.s) # select all if global test is not rejected
        all.s <- aperm(all.s, c(1, 3, 2))
      }
      
      all.s.all <- abind(all.s.all, all.s, along =  2) # store to array
    }
    
    pis <- sort(unique(c(as.matrix(pv.all)))) # implicit parameter for ROC
    pis <- pis[pis < 1] 
    
    r1 <- sapply(pis, function(pi) mean(pv.all[unstab,] <= pi)) # FPR
    r2 <- sapply(pis, function(pi) mean(pv.all[stab,] <= pi)) # TPR
    
    if (t == 1) {
      plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t],
           xlab = "False positive rate", ylab = "True positive rate", cex.lab = exp.text, lwd = exp.lines)
    } else {
      points(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
    }
    pi0 <- 0.01 # consider default value
    points(mean(pv.all[unstab,] <= pi0), mean(pv.all[stab,] <= pi0), pch = 4, col = cols[t], cex = exp.points)
    points(mean(all.s.all[,unstab,]), mean(all.s.all[,stab,]), pch = pchs[t], col = cols[t], cex = exp.points)
  }
  abline(0,1, col = "gray", lty= 2)
  labels <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns), "$')", sep = "", collapse = ","), ")")))
  legend('bottomright', legend = labels, col = cols, lty = ltys, pch = pchs, cex = exp.text, pt.cex = 1, lwd = exp.lines)
  
  # plotting parameters  
  cols <- 1:4
  ltys <- 1:4
  
  t <- 0
  ns.p <- numeric(0)
  for (file in flz){
    # loop through sample sizes
    t <- t + 1
    load(paste(folder, "/", file, sep = ""))
    # matrix and array to combine all information
    pv.corr.unstab <- matrix(NA, ncol = nsim, nrow = 0)
    pv.unstab <- array(NA, dim = c(nsim, B, 0))
    for (s in 1:ncol(all.comb)){
      # loop through covariate settings
      mes <- all.comb[, s]
      if (length(which.stab(mes)) < length(mes)){
        # find p-values for not well-specified settings
        pv <- simulation[[paste(mes, collapse = "")]]$pval
        pv.corr <- simulation[[paste(mes, collapse = "")]]$pval.corr
        pv.corr.unstab <- rbind(pv.corr.unstab, pv.corr) # store to matrix
        pv.unstab <- abind(pv.unstab, pv, along = 3) # store to array
      }
    }
    if(max(pv.unstab) < 5e-2) next() # stop plotting if always significant
    ns.p <- c(ns.p, simulation$n) # add sample size
    np <- prod(dim(pv.unstab)) # total number of p-values
    npc <- prod(dim(pv.corr.unstab)) # total number of combined p-values
    if (t == 1){
      # plot ECDF for all p-values
      plot(c(sort(pv.unstab), 1), (1:(np + 1))/(np + 1), type = "l", xlim = c(0, 1),
           col = cols[t], lty = ltys[t], xlab = "p", ylab = "Fn(p)",
           cex.lab = exp.text, lwd = exp.lines)
    } else {
      lines(c(sort(pv.unstab), 1), (1:(np + 1))/(np + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
    }
    # add ECDF for combined p-values
    lines(c(sort(pv.corr.unstab), 1), (1:(npc + 1))/(npc + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
    points(c(sort(pv.corr.unstab), 1), (1:(npc + 1))/(npc + 1), col = alpha(cols[t], 0.2))
  }
  labels.sub <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns.p), "$')", sep = "", collapse = ","), ")")))
  legend('bottomright', legend = labels.sub, col = cols, lty = ltys, cex = exp.text, pt.cex = 1, lwd = exp.lines)
  
  t <- 0
  ns.p <- numeric(0)
  for (file in flz){
    # loop through sample sizes
    t <- t + 1
    load(paste(folder, "/", file, sep = ""))
    # matrix and array to combine all information
    pv.corr.stab <- matrix(NA, ncol = nsim, nrow = 0)
    pv.stab <- array(NA, dim = c(nsim, B, 0))
    for (s in 1:ncol(all.comb)){
      # loop through covariate settings
      mes <- all.comb[, s]
      if (length(which.stab(mes)) >= length(mes)){
        # find p-values for well-specified settings
        pv <- simulation[[paste(mes, collapse = "")]]$pval
        pv.corr <- simulation[[paste(mes, collapse = "")]]$pval.corr
        pv.corr.stab <- rbind(pv.corr.stab, pv.corr) # store to matrix
        pv.stab <- abind(pv.stab, pv, along = 3) # store to array
      }
    }

    ns.p <- c(ns.p, simulation$n) # add sample size
    np <- prod(dim(pv.stab)) # total number of p-values
    npc <- prod(dim(pv.corr.stab)) # total number of combined p-values
    if (t == 1){
      # plot ECDF for all p-values
      plot(c(sort(pv.stab), 1), (1:(np + 1))/(np + 1), type = "l", xlim = c(0, 1),
           col = cols[t], lty = ltys[t], xlab = "p", ylab = "Fn(p)", cex.lab = exp.text, lwd = exp.lines)
    } else {
      lines(c(sort(pv.stab), 1), (1:(np + 1))/(np + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
    }
    # add ECDF for combined p-values
    lines(c(sort(pv.corr.stab), 1), (1:(npc + 1))/(npc + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
    points(c(sort(pv.corr.stab), 1), (1:(npc + 1))/(npc + 1), col = alpha(cols[t], 0.2))
  }
  labels.sub <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns.p), "$')", sep = "", collapse = ","), ")")))
  legend('bottomright', legend = labels.sub, col = cols, lty = ltys, cex = exp.text, pt.cex = 1, lwd = exp.lines)
}

figures6 <- function(folder){
  flz <- list.files(folder)
  flz <- flz[grep("results", flz)]
  
  
  # plotting parameters
  par(mfrow = c(1, 3))
  exp.text <- 1.5
  exp.points <- 1.5
  exp.lines <- 1.5
  cols <- ltys <-  1:4
  pchs <- c(0:2, 6)
  
  
  load(paste(folder, "/", flz[1], sep = ""))
  # simulation runs per sample size
  B <- dim(simulation$steps.out)[1] # number of splits used
  p <- dim(simulation$steps.out)[2]
  stab <- c(2)
  unstab <- (1:p)[-stab]
  t <- 0
  ns.p <- numeric(0)
  for (file in flz){
    # loop through sample sizes
    t <- t + 1
    load(paste(folder, "/", file, sep = ""))
    all.s.0 <- simulation$steps.out
    pi0 <- 1e-2
    
    if (get.pval){
      pval.lim <- 0.05
      pval <- simulation$pval.corr
      glob <- which(pval > pval.lim)
      glob.s <- t(simulation$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
    }
    
    pv <- apply(1 * is.na(all.s.0), 3, fisher.split)
    if(get.pval) pv[, glob] <- 0
    pis <- sort(unique(c(pv)))
    pis <- pis[pis < 1]
    
    r1 <- sapply(pis, function(pi) mean(pv[unstab,] <= pi))
    r2 <- sapply(pis, function(pi) mean(pv[stab,] <= pi))
    if (t == 1) {
      plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t],
           xlab = "False positive rate", ylab = "True positive rate", cex.lab = exp.text, lwd = exp.lines)
    } else {
      points(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
    }
    
    points(mean(pv[unstab,] <= pi0), mean(pv[stab,] <= pi0), pch = 4, col = cols[t], cex = exp.points)
    
    t1.mat <- is.na(all.s.0[,unstab,])
    p.mat <- is.na(all.s.0[,stab,])
    if (get.pval){
      if(length(dim(t1.mat)) > 2)
        t1.mat <- aperm(t1.mat, c(1, 3, 2))
      if(length(dim(p.mat)) > 2)
        p.mat <- aperm(p.mat, c(1, 3, 2))
      t1.mat <- pmax(t1.mat, glob.s)
      p.mat <- pmax(p.mat, glob.s)
    }
    
    
    points(mean(t1.mat), mean(p.mat), pch = pchs[t], col = cols[t], cex = exp.points)
    # points(mean(which.sel.tau[,unstab]), mean(which.sel.tau[,stab]), col = 4, pch = 4)
    abline(0,1, col = "gray", lty= 2)
  }
  labels <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns), "$')", sep = "", collapse = ","), ")")))
  legend('bottomright', legend = labels, col = cols, lty = ltys, pch = pchs, cex = exp.text, pt.cex = 1, lwd = exp.lines)
  
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
  
  s <- 0
  rdifs <- list()
  for (file in flz){
    if(grepl("results", file)){
      s <- s + 1
      load(paste(folder, "/", file, sep = ""))
      rdifs[[s]] <- simulation$all[,3] / ns[s]
    }
  }
  rdif <- matrix(unlist(rdifs), ncol = nf)
  boxplot(rdif, log = "y", names = labels, ylab = "Average misposition", cex.lab = exp.text, cex.axis =exp.text, yaxt ="n")
  axis(side = 2)
}


