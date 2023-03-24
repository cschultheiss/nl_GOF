d <- "/usr/local64.sfs/app/R/R_local/library"
if(all(d != .libPaths())) .libPaths(c(.libPaths(), d))

require(tictoc)
require(doRNG)
require(doSNOW)
require(git2r)
require(FOCI)
require(mgcv)
require(sfsmisc)

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
Ex <- function(x2) -(dnorm(up(x2)) - dnorm(lo(x2)))/(pnorm(up(x2)) - pnorm(lo(x2)))
Vx <- function(x2) 1 - (up(x2) * dnorm(up(x2)) - lo(x2) * dnorm(lo(x2)))/(pnorm(up(x2)) - pnorm(lo(x2))) -
  ((dnorm(up(x2)) - dnorm(lo(x2)))/(pnorm(up(x2)) - pnorm(lo(x2))))^2

nsim <- 200
n.vec <- 10^(2:5)
n.split <- 25
p <- 2
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
               .packages = c("mgcv"), .options.snow = opts) %dorng%{
                 
                 if(all(d != .libPaths())) .libPaths(c(.libPaths(), d))
                  library(FOCI)
                 
                 x0 <- rnorm(n)
                 x1 <- a * pot(x0 , b)/sqrt(va(b)) + runif(n, -1, 1)
                 x2 <- x1 + rnorm(n)
                 y <- x0 + pot(x2, 1.5)
                 dat <- data.frame(y, x1, x2)
                 form <- wrapFormula(y ~., data = dat)
                 fi.all <- gam(form, data = dat)
                 sel.all <- (1:2) %in% 
                   foci(abs(fi.all$residuals), dat[,-1])$selectedVar$index
                 sel.all0 <- (1:2) %in% 
                   foci(abs(x0 - Ex(x1)), dat[,-1])$selectedVar$index
                 
                 mse <- mean((fi.all$residuals - (x0 - Ex(x1)))^2)
                 rcor <- cor(fi.all$residuals, x0 - Ex(x1), method = "spearman")
                 sel11 <- sel12 <- sel21 <- sel22 <- matrix(NA, n.split, dim(dat)[2] - 1)
                 steps11 <-  matrix(NA, nrow = n.split, ncol = dim(dat)[2] - 1)
                 colnames(steps11) <- colnames(dat)[-1]
                 steps12 <- steps21 <- steps22 <- steps11
                 for(i in 1:n.split){
                   ind <- sample(n, n/2)
                   fi1 <- gam(form, data = dat[ind,])
                   fi2 <- gam(form, data = dat[-ind,])
                   pred11 <- predict(fi1, newdata = dat[ind,])
                   pred12 <- predict(fi1, newdata = dat[-ind,])
                   pred21 <- predict(fi2, newdata = dat[ind,])
                   pred22 <- predict(fi2, newdata = dat[-ind,])
                   fo11 <- foci(abs(y[ind] - pred11), dat[ind, -1])
                   fo12 <- foci(abs(y[-ind] - pred12), dat[-ind, -1])
                   fo21 <- foci(abs(y[ind] - pred21), dat[ind, -1])
                   fo22 <- foci(abs(y[-ind] - pred22), dat[-ind, -1])
                   steps11[i, fo11$selectedVar$index] <- diff(c(0, fo11$stepT))
                   steps12[i, fo12$selectedVar$index] <- diff(c(0, fo12$stepT))
                   steps21[i, fo21$selectedVar$index] <- diff(c(0, fo21$stepT))
                   steps22[i, fo22$selectedVar$index] <- diff(c(0, fo22$stepT))
                   sel11[i,] <- c(fo11$selectedVar$index, rep(NA, dim(dat)[2] - 1 - length(fo11$selectedVar$index)))
                   sel12[i,] <- c(fo12$selectedVar$index, rep(NA, dim(dat)[2] - 1 - length(fo12$selectedVar$index)))
                   sel21[i,] <- c(fo21$selectedVar$index, rep(NA, dim(dat)[2] - 1 - length(fo21$selectedVar$index)))
                   sel22[i,] <- c(fo22$selectedVar$index, rep(NA, dim(dat)[2] - 1 - length(fo22$selectedVar$index)))
                 }
                 
                 out <- list()
                 out$values <- c(mse, rcor, sel.all0, sel.all)
                 out$steps.in <- rbind(steps11, steps22)
                 out$steps.out <- rbind(steps21, steps12)
                 out$sel.in <- rbind(sel11, sel22)
                 out$sel.out <- rbind(sel21, sel12)
                 
                 out
               }
  toc()
  stopCluster(cl)
  
  # store output list to matrix
  res.val <- matrix(unlist(res[, "values"]), byrow = TRUE, nrow = nsim)
  colnames(res.val) <- c("mse", "rcor",
                         paste(rep(c("all0", "all"), each = p), rep(colnames(res[1,"steps.in"][[1]]), 2), sep ="."))
  
  res.steps.in <- array(unlist(res[,"steps.in"]), dim = c(2 * n.split, p, nsim), dimnames = list(NULL,
                                                                                 colnames(res[1,"steps.in"][[1]]), NULL))
  res.steps.out <- array(unlist(res[,"steps.out"]), dim = c(2 * n.split, p, nsim), dimnames = list(NULL,
                                                                                                 colnames(res[1,"steps.out"][[1]]), NULL))
  res.sel.in <- array(unlist(res[,"sel.in"]), dim = c(2 * n.split, p, nsim))
  res.sel.out <- array(unlist(res[,"sel.out"]), dim = c(2 * n.split, p, nsim))
  
  # store output quantities, sample size, random seed, commit
  simulation <- list(all = res.val, steps.in = res.steps.in, steps.out = res.steps.out, sel.in = res.sel.in,
                     sel.out = res.sel.out, n = n, r.seed = attr(res, "rng"), "commit" = commit)
  # create unique filename based on sample size and time
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  # save the file to the folder
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  print(apply(res.val, 2, mean))
  for (j in 1:p){
    print(paste(j, ": ", sum(res.sel.in == j, na.rm = TRUE), sep = ""))
    print(paste(j, ": ", sum(res.sel.out == j, na.rm = TRUE), sep = ""))
  }
}
