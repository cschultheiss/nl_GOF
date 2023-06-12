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


nsim <- 200
n.vec <- 10^(2:5)

n.split <- 25
p.tot <- 5
p.mes <- 3
p.an <- 5

stab <- function(sel){
  all.comb <- combn(1:p.tot, p.mes)
  sel.pos <- which(duplicated(rbind(sel, t(combn(1:5,3)))))- 1
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


RNGkind("L'Ecuyer-CMRG")
# make it reproducible
set.seed(42)
mes <- replicate(nsim, sort(sample.int(p.tot, p.mes)))
disti.pool <- c(rnorm, function(n) runif(n, -sqrt(3), sqrt(3)), 
                function(n) rexp(n) * sample(c(-1, 1), n, replace = TRUE) / sqrt(2))
disti <- replicate(nsim, sample(rep(disti.pool, each = ceiling((p.tot + 1) / length(disti.pool))), p.tot + 1))
pots <- replicate(nsim, cbind(runif(p.an, 0.5, 0.8), runif(p.an, 1.2, 1.5))[cbind(1:p.an,sample(1:2,p.an, TRUE))])

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
                 
                 eps <- matrix(NA, n, p.tot + 1)
                 for (j in 1:(p.tot + 1)) {
                   eps[,j] <- disti[j, gu][[1]](n)
                 }
                 x1 <- eps[,1]
                 x2 <- scale(pot(x1, pots[1, gu])) + 0.5 * eps[,2]
                 x3 <- scale(pot(x2, pots[2, gu])) + 0.5 * eps[,3]
                 y <- scale(pot(x1, pots[3, gu])) + scale(pot(x3, pots[4, gu])) + 0.5 * eps[,6]
                 x4 <- scale(pot(y, pots[5, gu])) + 0.5 * eps[,4]
                 x5 <- eps[,5]
                 x.all <- data.frame(x1, x2, x3, x4, x5)
                 x <- x.all[,mes[,gu]]
                 dat <- data.frame(y, x)

                 
                 out <- multi.spec(dat, B = n.split, return.predictor = TRUE, return.residual = FALSE,
                                   fitting = fitxg, predicting = predxg)
                 
                 
                 out
               }
  toc()
  stopCluster(cl)
  
  # store output list to matrix
  
  res.steps.out <- array(unlist(res[,"steps"]), dim = c(2 * n.split, p, nsim), dimnames = list(NULL,
                                                                                               colnames(res[1,"steps"][[1]]), NULL))
  res.sel.out <- array(unlist(res[,"sel"]), dim = c(2 * n.split, p, nsim))
  
  res.pval <- matrix(unlist(res[, "pval"]), byrow = TRUE, nrow = nsim)
  
  res.pval.corr <- unlist(res[,"pval.corr"])
  
  # store output quantities, sample size, random seed, commit
  simulation <- list(steps.out = res.steps.out, sel.out = res.sel.out,
                     pval = res.pval, pval.corr = res.pval.corr,
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
