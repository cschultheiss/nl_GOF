rm(list = ls(all = TRUE))
d <- "/usr/local64.sfs/app/R/R_local/library"
if(all(d != .libPaths())) .libPaths(c(.libPaths(), d))

require(tictoc)
require(doRNG)
require(doSNOW)
require(git2r)
require(FOCI)
require(mgcv)
require(sfsmisc)
require(dHSIC)
require(xgboost)

source("multi_spec.R")
source("fitpred.R")

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))

save <- TRUE
# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("results/", newdir, sep=""))
}

# update on simulation progress
progress <- function(n, tag) {
  mod <- 16
  if (n %% mod == 0 ) {
    cat(sprintf('tasks completed: %d; tag: %d\n', n, tag))
  }
  if (n %% mod == 0 ) {
    toc()
    tic()
  }
}

opts <- list(progress = progress)

pot <- function(x, b) sign(x)*abs(x)^b
va <- function(b) 2^(b)*gamma(b + 0.5)/sqrt(pi)
up <- function(x2) pot((x2 + 1)*sqrt(va(b))/a, 1/b)
lo <- function(x2) pot((x2 - 1)*sqrt(va(b))/a, 1/b)
fx2 <- function(x2) 0.5 *(pnorm(up(x2)) - pnorm(lo(x2)))
Ex1 <- function(x2) -(dnorm(up(x2)) - dnorm(lo(x2)))/(pnorm(up(x2)) - pnorm(lo(x2)))
# Vx <- function(x2) 1 - (up(x2) * dnorm(up(x2)) - lo(x2) * dnorm(lo(x2)))/(pnorm(up(x2)) - pnorm(lo(x2))) -
#   ((dnorm(up(x2)) - dnorm(lo(x2)))/(pnorm(up(x2)) - pnorm(lo(x2))))^2

Ex <- function(x1, x2, x6, x7) 0.5 * (x1^2 + x2^2 + 2) + 2 * (Ex1(x6) + pot(x7, 1.5))

nsim <- 200
n.vec <- 10^(2:5)
n.split <- 25
p <- 4
b <- 1.5
a <- sqrt(1/3)


RNGkind("L'Ecuyer-CMRG")
# make it reproducible
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0


for (n in n.vec) {
  seed.n <- seed.n + 1
  # use different known seed for different n => can ommit lower sample size if wished
  set.seed(seed.vec[seed.n])
  
  # initiliaze parralelisation
  cl<-makeSOCKcluster(16)
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("mgcv", "sfsmisc", "xgboost"), .options.snow = opts) %dorng%{
                 
                 if(all(d != .libPaths())) .libPaths(c(.libPaths(), d))
                   library(FOCI)
                   library(dHSIC)
                
                
                 x0 <- rnorm(n)
                 x1 <- sqrt(0.5) * (x0 + rnorm(n))
                 x2 <- sqrt(0.5) * (x0 + rnorm(n))
                 x3 <- sqrt(0.5) * (x1 + rnorm(n))
                 x4 <- sqrt(0.5) * (x2 + rnorm(n))
                 x5 <- rnorm(n)
                 x6 <- a * pot(x5 , b)/sqrt(va(b)) + runif(n, -1, 1)
                 x7 <- x6 + rnorm(n)
                 y <- x3^2 + x4^2 + 2 * (x5 + pot(x7, 1.5)) + rnorm(n, sd = 1)
                 
                 Eyx <- x3^2 + x4^2 + 2 * (x5 + pot(x7, 1.5))
                 
                 dat <- data.frame(y, x3, x4, x5, x7)
                 
                 eps0 <- y-Eyx
                 
                 steps.all <- steps.all0 <- rep(NA, p)
                 fi.all <- gam(wrapFormula(y ~ ., data = dat), data = dat)
                 foci.all <- foci(fi.all$residuals, dat[,-1])
                 steps.all[foci.all$selectedVar$index] <- diff(c(0, foci.all$stepT))
                 sel.all <- sapply(1:p, function(j) {
                   if (j %in% foci.all$selectedVar$index)
                     which(foci.all$selectedVar$index == j)
                   else
                     NA})
                 foci.all0 <- foci(eps0, dat[,-1])
                 steps.all0[foci.all0$selectedVar$index] <- diff(c(0, foci.all0$stepT))
                 sel.all0 <- sapply(1:p, function(j) {
                   if (j %in% foci.all0$selectedVar$index)
                     which(foci.all0$selectedVar$index == j)
                   else
                     NA})

                 
                 div <- which(abs(fi.all$residuals) > (max(abs(fi.all$residuals)) - 1e-5))
                 if(length(div) > 1){
                   mse <- mean((fi.all$residuals - eps0)[-div]^2)
                 } else {
                   mse <- mean((fi.all$residuals - eps0)^2)
                 }
                 
                 rcor <- cor(fi.all$residuals, eps0, method = "spearman")
                 rdif <- mean(abs(rank(fi.all$residuals) - rank(eps0)))
                 
                 out <- multi.spec(dat, B = n.split, return.predictor = FALSE, return.residual = TRUE)
                 out$co0 <- cor(out$residual, eps0, method = "spearman")
                 out$co <- cor(out$residual, method = "spearman")
                 diag(out$co) <- NA
                 out$co <- apply(out$co, 1, median, na.rm = TRUE)
                 out$residual = NULL

                 out$values <- c(mse, rcor, rdif, length(div), sel.all0, sel.all, steps.all0, steps.all)

                 out
               }
  toc()
  stopCluster(cl)
  
  # store output list to matrix
  res.val <- matrix(unlist(res[, "values"]), byrow = TRUE, nrow = nsim)
  colnames(res.val) <- c("mse", "rcor", "rdif", "0div",
                         paste(rep(c("pos0", "pos"), each = p), rep(colnames(res[1,"steps"][[1]]), 2), sep ="."),
                         paste(rep(c("all0", "all"), each = p), rep(colnames(res[1,"steps"][[1]]), 2), sep ="."))
  
  res.steps.out <- array(unlist(res[,"steps"]), dim = c(2 * n.split, p, nsim), dimnames = list(NULL,
                                                                                                 colnames(res[1,"steps"][[1]]), NULL))
  res.sel.out <- array(unlist(res[,"sel"]), dim = c(2 * n.split, p, nsim))
  
  res.pval <- matrix(unlist(res[, "pval"]), byrow = TRUE, nrow = nsim)
  
  res.co0 <- matrix(unlist(res[, "co0"]), byrow = TRUE, nrow = nsim)
  
  res.co <- matrix(unlist(res[, "co"]), byrow = TRUE, nrow = nsim)
  
  res.pval.corr <- unlist(res[,"pval.corr"])
  
  # store output quantities, sample size, random seed, commit
  simulation <- list(all = res.val, steps.out = res.steps.out, sel.out = res.sel.out,
                     pval = res.pval, pval.corr = res.pval.corr, co0 = res.co0, co = res.co,
                     n = n, r.seed = attr(res, "rng"), "commit" = commit)
  # create unique filename based on sample size and time
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  # save the file to the folder
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  print(median(res.pval.corr))
  print(apply(res.val[,1:4], 2, mean))
  print(apply(!is.na(res.val[,(4 + 1) : (4 + 2 * p)]), 2, mean))
  for (j in 1:p){
    # print(paste(j, ": ", sum(res.sel.in == j, na.rm = TRUE), sep = ""))
    print(paste(j, ": ", sum(res.sel.out == j, na.rm = TRUE), sep = ""))
  }
}
