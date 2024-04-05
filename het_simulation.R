rm(list = ls(all = TRUE))
het_simulation <- function(nsim = 200, n.vec = 10^(2:5), extra.packages = NULL, extra.path = NULL){
  # when called executes the simulation for Section 6.2
  # Input
  # nsim (integer): number of simulation runs per sample size
  # n.vec (integer vector): sample sizes to be considered
  # extra.packages (character, vector): packages not defined in usual location
  # extra.path (character): path to find these packages
  # Output (string): save location of the results folder
  if(!is.null(extra.path)){
    d <- extra.path
    if(all(d != .libPaths())) .libPaths(c(.libPaths(), d))
  }
  
  export.packages <- c("mgcv", "sfsmisc", "xgboost", "FOCI", "dHSIC")
  export.packages <- setdiff(export.packages, extra.packages)
  
  require(tictoc)
  require(doRNG)
  require(doSNOW)
  require(git2r)
  require(FOCI)
  require(mgcv)
  require(sfsmisc)
  require(dHSIC)
  require(xgboost)
  
  source("multi_spec.R", local = TRUE)
  source("fitpred.R", local = TRUE)
  
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
  
  
  fh <- function(h) dnorm(h) # distribution of H
  fx1 <- function(x) dnorm(x, sd = sqrt(0.5)) # distribution of xi_1
  f1 <- function(h) sin(2 * h) # effect from H to X_1
  fhx1 <- function(h, x) fh(h) * fx1(x - f1(h)) # joint distribution of H and X_1
  # marginal distribution of X_1
  norm <- Vectorize(function(x) integrate(function(h) fhx1(h, x), 2 * min(h), 2 * max(h))$value)
  # first conditional moment E[H|X_1] up to normalization
  m1 <- Vectorize(function(x) integrate(function(h) h * fhx1(h, x), 2 * min(h), 2 * max(h))$value)
  # first conditional moment E[H|X_2] up to normalization
  m2 <- Vectorize(function(x) integrate(function(h) h^2 * fhx1(h, x), 2 * min(h), 2 * max(h))$value)
  
  RNGkind("L'Ecuyer-CMRG")
  # make it reproducible
  # sample a large number to find range
  set.seed(42)
  n <- 1e7
  h <- rnorm(n)
  eps <- rnorm(n, sd = sqrt(0.5))
  x1 <- f1(h) + eps 
  # restrict to smaller grid
  x.grid <- seq(quantile(x1, 0), quantile(x1, 1), length.out = 1e3)
  # remove the large vectors
  rm(x1, eps)
  h <- c(min(h), max(h))
  norm.grid <- norm(x.grid)
  # enforce unimodality => more stable estimation
  x.max <- which.max(norm.grid)
  lower <- which(norm.grid[1:x.max] < cummax(norm.grid[1:x.max]))
  higher <- which(norm.grid[x.max:1e3] > cummin(norm.grid[x.max:1e3]))
  omit <- c(lower, higher)
  if(length(omit) > 0){
    x.grid <- x.grid[-omit]
    norm.grid <- norm.grid[-omit]
  }
  # fit splines for faster evaluation downstream
  ss1 <- smooth.spline(x.grid, m1(x.grid)/norm.grid, cv = TRUE)
  ss2 <- smooth.spline(x.grid, m2(x.grid)/norm.grid, cv = TRUE)
  
  n.split <- 25 # number of splits to assess well-specifaction
  p <- 2 # number of covariates
  
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
                 .packages = export.packages, .options.snow = opts) %dorng%{
                   
                   # some packages must be load on the workers if not defined in default environment
                   if(!is.null(extra.packages)){
                     if(all(d != .libPaths())) .libPaths(c(.libPaths(), d))
                     lapply(extra.packages, require, character.only = TRUE)
                   }
                   
                   h <- rnorm(n) # simulate h
                   x1 <- f1(h) + rnorm(n, sd = sqrt(0.5)) # simulate x1
                   x2 <- (x1 + rnorm(n, sd = sqrt(0.5)))*sqrt(2/3) # simulate x2
                   
                   y <- (abs(x2^2) + 2)/(abs(x2^2) + 1) * (h) # simulate 1
                   
                   m1 <- (abs(x2^2) + 2)/(abs(x2^2) + 1) * (predict(ss1, x1)$y) # calculate first moment
                   m2 <- ((abs(x2^2) + 2)/(abs(x2^2) + 1))^2 *
                     (predict(ss2, x1)$y) # calculate second moment
                   eps0 <- (y - m1)/sqrt(m2 - m1^2) # calculate true residual
                   
                   dat <- data.frame(y, x1, x2)
                   steps.all <- steps.all0 <- rep(NA, p)
                   fi.all <- allfitxg.het(dat) # fit on all data
                   foci.all <- foci(fi.all$residuals, dat[,-1]) # apply FOCI once
                   steps.all[foci.all$selectedVar$index] <- diff(c(0, foci.all$stepT))
                   sel.all <- sapply(1:p, function(j) {
                     if (j %in% foci.all$selectedVar$index)
                       which(foci.all$selectedVar$index == j)
                     else
                       NA})
                   foci.all0 <- foci(eps0, dat[,-1]) # apply FOCI with true residual
                   steps.all0[foci.all0$selectedVar$index] <- diff(c(0, foci.all0$stepT))
                   sel.all0 <- sapply(1:p, function(j) {
                     if (j %in% foci.all0$selectedVar$index)
                       which(foci.all0$selectedVar$index == j)
                     else
                       NA})
                   
                   # find extreme residua, i.e., likely negative normalization
                   div <- which(abs(fi.all$residuals) > (max(abs(fi.all$residuals)) - 1e-5))
                   if(length(div) > 0){
                     mse <- mean((fi.all$residuals - eps0)[-div]^2)
                   } else {
                     mse <- mean((fi.all$residuals - eps0)^2)
                   }
                   
                   rcor <- cor(fi.all$residuals, eps0, method = "spearman") # rank correlation
                   rdif <- mean(abs(rank(fi.all$residuals) - rank(eps0))) # average misposition
                   
                   out <- multi.spec(dat, B = n.split, return.predictor = FALSE, return.residual = FALSE,
                                     fitting = fitxg.het, predicting = predxg.het, norming = normxg.het,
                                     trafo = function(x) x) # assess well-specification
                   
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
    # print some intermediate results
    print(median(res.pval.corr))
    print(apply(res.val[,1:4], 2, mean))
    print(apply(!is.na(res.val[,(4 + 1) : (4 + 2 * p)]), 2, mean))
    for (j in 1:p){
      print(paste(j, ": ", sum(res.sel.out == j, na.rm = TRUE), sep = ""))
    }
  }
  return(paste("results/", newdir, sep = ""))
}

het_simulation(extra.packages = c("FOCI", "dHSIC"), extra.path = "/usr/local64.sfs/app/R/R_local/library")
