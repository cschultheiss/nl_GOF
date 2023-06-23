require(latex2exp)
require(abind)
require(scales)
source("split.R")

folder <- "results/rand-setup-gam"
savefolder <- "Figures/rand"
flz <- list.files(folder)
flz <- flz[grep("results", flz)]
p.tot <- 5
p.mes <- 3
load(paste(folder, "/", flz[1], sep = ""))
nsim <- dim(simulation[[4]]$pval)[1]
B <- dim(simulation[[4]]$pval)[2]

all.comb <- combn(1:p.tot, p.mes)

which.stab <- function(sel){
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

stab.col <- rep(FALSE, prod(dim(all.comb)))

for (s in 1:ncol(all.comb)){
  mes <- all.comb[, s]
  stab.col[(s - 1) * p.mes + which(mes %in% which.stab(mes))] <- TRUE
  if(length(which.stab(mes)) == length(mes))
    stab.col[(s - 1) * p.mes + 1:3] <- NA
    
}
stab <- which(stab.col)
unstab <- which(!stab.col)

get.pval <- TRUE

# png(paste(savefolder, "/ROC+ECDF.png", sep = ""), width = 3600,
# height = 1200, res = 300)
par(mfrow = c(1, 3))
exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
t <- 0
ns <- numeric(length(flz))
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  ns[t] <- simulation$n
  pv.all <- matrix(NA, ncol = nsim, nrow = 0)
  all.s.all <- array(NA, dim = c(B, 0, nsim))
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    all.s <- simulation[[paste(mes, collapse = "")]]$steps.out
    if (get.pval){
      pval.lim <- 0.05
      glob <- which(simulation[[paste(mes, collapse = "")]]$pval.corr > pval.lim)
      glob.s <- t(simulation[[paste(mes, collapse = "")]]$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
    }
    pv <- apply(1 * is.na(all.s), 3, fisher.split)
    if(get.pval) pv[, glob] <- 0
    pv.all <- rbind(pv.all, pv)
    
    all.s <- is.na(all.s)
    if (get.pval){
      all.s <- aperm(all.s, c(1, 3, 2))
      all.s <- pmax(all.s, glob.s)
      all.s <- aperm(all.s, c(1, 3, 2))
    }
    
    all.s.all <- abind(all.s.all, all.s, along =  2)
  }
  
  pis <- sort(unique(c(as.matrix(pv.all))))
  pis <- pis[pis < 1]
  
  r1 <- sapply(pis, function(pi) mean(pv.all[unstab,] <= pi))
  r2 <- sapply(pis, function(pi) mean(pv.all[stab,] <= pi))
  
  if (t == 1) {
    plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t],
         xlab = "False positive rate", ylab = "True positive rate", cex.lab = exp.text, lwd = exp.lines)
  } else {
    points(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
  pi0 <- 0.01
  points(mean(pv.all[unstab,] <= pi0), mean(pv.all[stab,] <= pi0), pch = 4, col = cols[t], cex = exp.points)
  points(mean(all.s.all[,unstab,]), mean(all.s.all[,stab,]), pch = pchs[t], col = cols[t], cex = exp.points)
}
abline(0,1, col = "gray", lty= 2)
labels <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns), "$')", sep = "", collapse = ","), ")")))
legend('bottomright', legend = labels, col = cols, lty = ltys, pch = pchs, cex = exp.text, pt.cex = 1, lwd = exp.lines)

cols <- 1:4
ltys <- 1:4
t <- 0
ns.p <- numeric(0)
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  pv.corr.unstab <- matrix(NA, ncol = nsim, nrow = 0)
  pv.unstab <- array(NA, dim = c(nsim, B, 0))
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    if (length(which.stab(mes)) < length(mes)){
    pv <- simulation[[paste(mes, collapse = "")]]$pval
    pv.corr <- simulation[[paste(mes, collapse = "")]]$pval.corr
    pv.corr.unstab <- rbind(pv.corr.unstab, pv.corr)
    pv.unstab <- abind(pv.unstab, pv, along = 3)
    }
  }
  if(max(pv.unstab) < 5e-2) next()
  ns.p <- c(ns.p, simulation$n)
  np <- prod(dim(pv.unstab))
  npc <- prod(dim(pv.corr.unstab))
  if (t == 1){
    plot(c(sort(pv.unstab), 1), (1:(np + 1))/(np + 1), type = "l", xlim = c(0, 1),
         col = cols[t], lty = ltys[t], xlab = "p", ylab = "Fn(p)",
         cex.lab = exp.text, lwd = exp.lines)
  } else {
    lines(c(sort(pv.unstab), 1), (1:(np + 1))/(np + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
  
  lines(c(sort(pv.corr.unstab), 1), (1:(npc + 1))/(npc + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
  points(c(sort(pv.corr.unstab), 1), (1:(npc + 1))/(npc + 1), col = alpha(cols[t], 0.2))
}
labels.sub <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns.p), "$')", sep = "", collapse = ","), ")")))
legend('bottomright', legend = labels.sub, col = cols, lty = ltys, cex = exp.text, pt.cex = 1, lwd = exp.lines)

t <- 0
ns.p <- numeric(0)
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  pv.corr.stab <- matrix(NA, ncol = nsim, nrow = 0)
  pv.stab <- array(NA, dim = c(nsim, B, 0))
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    if (length(which.stab(mes)) >= length(mes)){
      pv <- simulation[[paste(mes, collapse = "")]]$pval
      pv.corr <- simulation[[paste(mes, collapse = "")]]$pval.corr
      pv.corr.stab <- rbind(pv.corr.stab, pv.corr)
      pv.stab <- abind(pv.stab, pv, along = 3)
    }
  }
  # print(c(ks.test(pv.stab, punif, alternative = "greater")$p.value,
  #         ks.test(pv.corr.stab, punif, alternative = "greater")$p.value))
  if(max(pv.stab) < 5e-2) next()
  ns.p <- c(ns.p, simulation$n)
  np <- prod(dim(pv.stab))
  npc <- prod(dim(pv.corr.stab))
  if (t == 1){
    plot(c(sort(pv.stab), 1), (1:(np + 1))/(np + 1), type = "l", xlim = c(0, 1),
         col = cols[t], lty = ltys[t], xlab = "p", ylab = "Fn(p)", cex.lab = exp.text, lwd = exp.lines)
  } else {
    lines(c(sort(pv.stab), 1), (1:(np + 1))/(np + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
  
  lines(c(sort(pv.corr.stab), 1), (1:(npc + 1))/(npc + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
  points(c(sort(pv.corr.stab), 1), (1:(npc + 1))/(npc + 1), col = alpha(cols[t], 0.2))
}
labels.sub <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns.p), "$')", sep = "", collapse = ","), ")")))
legend('bottomright', legend = labels.sub, col = cols, lty = ltys, cex = exp.text, pt.cex = 1, lwd = exp.lines)
# dev.off()

par(mfrow = c(2,4))
exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
t <- 0
ns <- numeric(length(flz))
all.comb.res <- all.comb[, c(1,5,3,6,8:10,2,4,7)]
for (s in 1:ncol(all.comb)){
  t <- 0
  mes <- all.comb.res[, s]
  stab.col <- which.stab(mes)
  stab.s <- which(mes %in% stab.col)
  unstab.s <- which(!(mes %in% stab.col))
  if(length(stab.col) == p.mes) next()
  for (file in flz){
    load(paste(folder, "/", file, sep = ""))
    t <- t + 1
    ns[t] <- simulation$n
    all.s <- simulation[[paste(mes, collapse = "")]]$steps.out
    if (get.pval){
      pval.lim <- 0.05
      glob <- which(simulation[[paste(mes, collapse = "")]]$pval.corr > pval.lim)
      glob.s <- t(simulation[[paste(mes, collapse = "")]]$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
    }
    pv <- apply(1 * is.na(all.s), 3, fisher.split)
    if(get.pval) pv[, glob] <- 0
    
    all.s <- is.na(all.s)
    if (get.pval){
      all.s <- aperm(all.s, c(1, 3, 2))
      all.s <- pmax(all.s, glob.s)
      all.s <- aperm(all.s, c(1, 3, 2))
    }
    
    pis <- sort(unique(c(as.matrix(pv))))
    pis <- pis[pis < 1]
    
    r1 <- sapply(pis, function(pi) mean(pv[unstab.s,] <= pi))
    r2 <- sapply(pis, function(pi) mean(pv[stab.s,] <= pi))
    
    if(length(stab.s) > 0){
      if (t == 1) {
        plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t],
             xlab = "False positive rate", ylab = "True positive rate", main = paste(mes, collapse = ", "),
             cex.lab = exp.text, lwd = exp.lines)
      } else {
        points(r1, r2, type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
      }
      pi0 <- 0.01
      points(mean(pv[unstab.s,] <= pi0), mean(pv[stab.s,] <= pi0), pch = 4, col = cols[t], cex = exp.points)
      points(mean(all.s[,unstab.s,]), mean(all.s[,stab.s,]), pch = pchs[t], col = cols[t], cex = exp.points)
    abline(0,1, col = "gray", lty= 2)
    } else {
      if (t == 1) {
        plot(r1, pis, xlim = c(0,1), ylim = c(0,0.5), type = "l", col = cols[t], lty = ltys[t],
             xlab = "False positive rate", ylab = "Threshold", main = paste(mes, collapse = ", "),
             cex.lab = exp.text, lwd = exp.lines)
      } else {
        points(r1, pis, type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
      }
      pi0 <- 0.01
      points(mean(pv[unstab.s,] <= pi0), pi0, pch = 4, col = cols[t], cex = exp.points)
      points(mean(all.s[,unstab.s,]), max(pis), pch = pchs[t], col = cols[t], cex = exp.points)
    }
  }
  if (s== 4){
    labels <- eval(parse(text = paste("c(", paste("TeX('$n=10^", log10(ns), "$')", sep = "", collapse = ","), ")")))
    legend('bottomright', legend = labels, col = cols, lty = ltys, pch = pchs, cex = exp.text, pt.cex = 1, lwd = exp.lines)
  }  
}


exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
get.pval <- FALSE
t <- 0
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  pv0.all <- matrix(NA, ncol = nsim, nrow = 0)
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    if(length(which.stab(mes)) < length(mes)) next()
    pv0 <- simulation[[paste(mes, collapse = "")]]$H0
    pv0.all <- rbind(pv0.all, pv0)
  } 
  print(ks.test(pv0.all, punif)$p.value)
  np <- prod(dim(pv0.all))
  if (t == 1){
    plot(c(sort(pv0.all), 1), (1:(np + 1))/(np + 1), type = "l", xlim = c(0, 1),
         col = cols[t], lty = ltys[t], xlab = "p", ylab = "Fn(p)", cex.lab = exp.text, lwd = exp.lines)
  } else {
    lines(c(sort(pv0.all), 1), (1:(np + 1))/(np + 1), col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
}


exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
t <- 0
ns <- numeric(length(flz))
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  ns[t] <- simulation$n
  all.s.all <- array(NA, dim = c(B, 0, nsim))
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    all.s <- simulation[[paste(mes, collapse = "")]]$steps.out
    all.s.all <- abind(all.s.all, all.s, along =  2)
  }
  ct <- apply(is.na(all.s.all), 2:3, sum)
  
  if (t == 1) {
    qqplot(ct[stab,], ct[unstab,], xlim = c(0, B), ylim = c(0, B),
           type = "l", col = cols[t], lty = ltys[t], cex.lab = exp.text, lwd = exp.lines,
           xlab = "Well-specified", ylab = "Not well-specified")
  } else {
    le <- nsim * length(stab) - 1
    lines(sort(ct[stab,]), quantile(ct[unstab,], (0:le)/le), type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
}



exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
get.pval <- FALSE
t <- 0
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  pv.all <- matrix(NA, ncol = nsim, nrow = 0)
  all.s.all <- array(NA, dim = c(B, 0, nsim))
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    if(length(which.stab(mes)) < length(mes)) next()
    all.s <- simulation[[paste(mes, collapse = "")]]$steps.out
    B <- dim(all.s)[1]
    if (get.pval){
      pval.lim <- 0.05
      glob <- which(simulation[[paste(mes, collapse = "")]]$pval.corr > pval.lim)
      glob.s <- t(simulation[[paste(mes, collapse = "")]]$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
    }
    pv <- apply(1 * is.na(all.s), 3, fisher.split)
    if(get.pval) pv[, glob] <- 0
    pv.all <- rbind(pv.all, pv)
    
    all.s <- is.na(all.s)
    if (get.pval){
      all.s <- aperm(all.s, c(1, 3, 2))
      all.s <- pmax(all.s, glob.s)
      all.s <- aperm(all.s, c(1, 3, 2))
    }
    
    all.s.all <- abind(all.s.all, all.s, along =  2)
  }
  
  pis <- sort(unique(c(as.matrix(pv.all))))
  pis <- pis[pis < 1]
  

  r2 <- sapply(pis, function(pi) mean(pv.all <= pi))
  
  if (t == 1) {
    plot(pis, r2, xlim = c(0,max(pis)), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t],
         xlab = "Threshold", ylab = "True positive rate", cex.lab = exp.text, lwd = exp.lines)
  } else {
    points(pis, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
  pi0 <- 0.01
  points(pi0, mean(pv.all <= pi0), pch = 4, col = cols[t], cex = exp.points)
  
  points(0.5, mean(all.s.all), pch = pchs[t], col = cols[t], cex = exp.points)
}