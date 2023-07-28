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
require(causaldata)
require(drf)

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
n.split <- 25
cc <- close_college[,c("educ", "exper", "lwage")]
p <- ncol(cc)
y.col <- 4

n.vec <- c(10^3, nrow(cc), 10^(4:5))

RNGkind("L'Ecuyer-CMRG")
# make it reproducible
set.seed(42)
seed.vec <- sample(1:10000, length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0

DRF <- drf(cc[,c("educ", "exper")], cc$lwage)
weights <- get_sample_weights(DRF)

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
                 
                 ind <- sample.int(nrow(cc), n, replace = TRUE)
                 dat <- data.frame(cc[ind, c("educ", "exper")])
                 dat$lwage <- sapply(ind, function(i) sample(cc$lwage, 1, prob = weights[i,]))
                 
                 eps0 <- (rchisq(n, df = 1) - 1) / sqrt(8) * sd(cc$lwage) * sign(dat$lwage - median(cc$lwage))
                 dat$extra <- dat$lwage + eps0
                 colnames(dat)[y.col] <- "y"
                 
                 steps.all <- steps.all0 <- rep(NA, p)
                 fi.all <- gam(wrapFormula(y ~., data = dat), data = dat)
                 foci.all <- foci(fi.all$residuals, dat[,-y.col])
                 steps.all[foci.all$selectedVar$index] <- diff(c(0, foci.all$stepT))
                 sel.all <- sapply(1:p, function(j) {
                   if (j %in% foci.all$selectedVar$index)
                     which(foci.all$selectedVar$index == j)
                   else
                     NA})
                 foci.all0 <- foci(eps0, dat[,-y.col])
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
                 
                 out <- multi.spec(dat, B = n.split, return.predictor = FALSE, return.residual = FALSE)
                 
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
  
  res.pval.corr <- unlist(res[,"pval.corr"])
  
  # store output quantities, sample size, random seed, commit
  simulation <- list(all = res.val, steps.out = res.steps.out, sel.out = res.sel.out,
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
