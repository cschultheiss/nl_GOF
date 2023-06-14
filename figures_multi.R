require(latex2exp)
source("split.R")

folder <- "results/14-Jun-2023 11.43"
flz <- list.files(folder)
flz <- flz[grep("results", flz)]
p.tot <- 5
p.mes <- 3

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
    c(2, 5)
  } else {
    error("did not find a match")
  }
}

stab.col <- rep(FALSE, prod(dim(all.comb)))

for (s in 1:ncol(all.comb)){
  mes <- all.comb[, s]
  stab.col[(s - 1) * p.mes + which(mes %in% which.stab(mes))] <- TRUE
}
stab <- (1:prod(dim(all.comb)))[stab.col]
unstab <- (1:prod(dim(all.comb)))[-stab]

get.pval <- TRUE

par(mfrow = c(1, 3))
exp.text <- 1.5
exp.points <- 1.5
exp.lines <- 1.5
cols <- ltys <-  1:4
pchs <- c(0:2, 6)
t <- 0
for (file in flz){
  t <- t + 1
  load(paste(folder, "/", file, sep = ""))
  pv.all <- matrix(NA, ncol = 100, nrow = 0)
  for (s in 1:ncol(all.comb)){
    mes <- all.comb[, s]
    all.s <- simulation[[paste(mes, collapse = "")]]$steps.out
    B <- dim(all.s)[1]
    if (get.pval){
      pval.lim <- 0.05
      glob <- which(simulation[[paste(mes, collapse = "")]]$pval.corr > pval.lim)
      glob.s <- t(simulation[[paste(mes, collapse = "")]]$pval[,c((B/2 + 1) : B, 1 : (B/2))]) > pval.lim
    }
    pv <- apply(1 * is.na(all.s), 3, fisher.split)
    if(get.pval) pv[, glob] <- 0
    print(length(glob))
    
    pv.all <- rbind(pv.all, pv)
  }
  
  pis <- sort(unique(c(as.matrix(pv.all))))
  
  r1 <- sapply(pis, function(pi) mean(pv.all[unstab,] <= pi))
  r2 <- sapply(pis, function(pi) mean(pv.all[stab,] <= pi))
  
  if (t == 1) {
    plot(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t],
         xlab = "False positive rate", ylab = "True positive rate", cex.lab = exp.text, lwd = exp.lines)
  } else {
    points(r1, r2, xlim = c(0,1), ylim = c(0,1), type = "l", col = cols[t], lty = ltys[t], lwd = exp.lines)
  }
}
