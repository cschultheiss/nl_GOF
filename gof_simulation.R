require(tictoc)
require(doRNG)
require(doSNOW)
require(git2r)
require(FOCI)
require(mgcv)

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
n.vec <- 10^(3:4)
n.split <- 10
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
               .packages = c("FOCI", "mgcv"), .options.snow = opts) %dorng%{
                 
                 x0 <- rnorm(n)
                 x1 <- a * pot(x0 , b)/sqrt(va(b)) + runif(n, -1, 1)
                 x2 <- x1 + rnorm(n)
                 y <- x0 + pot(x2, 1.5)
                 dat <- data.frame(y, x1, x2)
                 fi.all <- gam(y ~ s(x1) + s(x2), data = dat)
                 sel.all <- (1:2) %in% 
                   foci(abs(fi.all$residuals), dat[,-1])$selectedVar$index
                 sel.all0 <- (1:2) %in% 
                   foci(abs(x0 - Ex(x1)), dat[,-1])$selectedVar$index
                 sel <- matrix(NA, n.split, dim(dat)[2] - 1)
                 mse <- mean((fi.all$residuals - (x0 - Ex(x1)))^2)
                 rcor <- cor(fi.all$residuals, x0 - Ex(x1), method = "spearman")
                 for(i in 1:10){
                   ind <- sample(n, n/2)
                   fi <- gam(y ~ s(x1) + s(x2), data = dat[-ind,])
                   pred <- predict(fi, newdata = dat[ind,])
                   sel[i,] <- (1:2) %in% 
                     foci(abs(y[ind] - pred), dat[ind, -1])$selectedVar$index
                 }
                 
                 out <- list()
                 out$values <- c(mse, rcor, sel.all0, sel.all, colMeans(sel))     
                 out
               }
  toc()
  stopCluster(cl)
  
  # store output list to matrix
  res.val <- matrix(unlist(res[, "values"]), byrow = TRUE, nrow = nsim)
  colnames(res.val) <- c("mse", "rcor",
                         paste(rep(c("all0", "all", "split"), each = 2), rep(c("x1", "x2"), 3), sep ="."))
  
  # store output quantities, sample size, random seed, commit
  simulation <- list(results = res.val,
                     n = n, r.seed = attr(res, "rng"), "commit" = commit)
  # create unique filename based on sample size and time
  resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  # save the file to the folder
  if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  print(apply(res.val, 2, mean))

}
