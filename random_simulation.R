rm(list = ls(all = TRUE))
rand_simulation <- function(nsim = 100, n.vec = 10^(2:5), extra.packages = NULL, extra.path = NULL){
  # when called executes the simulation for Section 4
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
  
  n.split <- 25 # number of splits to assess well-specification 
  p.tot <- 5 # total number of covariates
  p.mes <- 3 # number of observed covariates
  p.an <- 5 # number of edges
  
  pot <- function(x, p) abs(x)^p * sign(x) # monotonic non-lineartiy
  
  all.comb <- combn(1:p.tot, p.mes) # possible observed subsets
  stab <- function(sel){
    # variables with well-specified effect per setting
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
  
  
  RNGkind("L'Ecuyer-CMRG")
  # make it reproducible
  set.seed(42)
  # possible noise distribution
  disti.pool <- c(rnorm, function(n) runif(n, -sqrt(3), sqrt(3)), 
                  function(n) rexp(n) * sample(c(-1, 1), n, replace = TRUE) / sqrt(2))
  # select distributions
  disti <- replicate(nsim, sample(rep(disti.pool, each = ceiling((p.tot + 1) / length(disti.pool))), p.tot + 1))
  potsl <- replicate(nsim, runif(p.an, 0.5, 0.8)) # sample monotonic effects
  potsh <- replicate(nsim, runif(p.an, 1.2, 1.5)) # sample nonmonotonic effects
  frac <- replicate(nsim, runif(p.an, 0.3, 0.7)) # sample fractions
  
  setup <- list(dist = disti, potsl = potsl, potsh = potsh, frac = frac) # save the sampled settings
  if (save) save(setup, file = paste("results/", newdir, "/setup.RData", sep = ""))
  
  set.seed(42) # make it reproducible
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
                   
                   eps <- matrix(NA, n, p.tot + 1)
                   # sample noise terms
                   for (j in 1:(p.tot + 1)) {
                     eps[,j] <- disti[j, gu][[1]](n)
                   }
                   # causal assignments
                   x1 <- eps[,1]
                   x2 <- scale(frac[1, gu] * pot(x1, potsl[1, gu]) 
                               + (1 - frac[1, gu]) * abs(x1)^potsh[1, gu]) + 0.5 * eps[,2]
                   x3 <- scale(frac[2, gu] * pot(x2, potsl[2, gu]) 
                               + (1 - frac[2, gu]) * abs(x2)^potsh[2, gu]) + 0.5 * eps[,3]
                   y <- scale(-frac[3, gu] * pot(x1, potsl[3, gu]) 
                              + (1 - frac[3, gu]) * abs(x1)^potsh[3, gu]) + 
                     scale(frac[4, gu] * pot(x3, potsl[4, gu]) 
                           + (1 - frac[4, gu]) * abs(x3)^potsh[4, gu]) + 0.5 * eps[,6]
                   x4 <- scale(frac[5, gu] * pot(y, potsl[5, gu]) 
                               + (1 - frac[5, gu]) * abs(y)^potsh[5, gu]) + 0.5 * eps[,4]
                   x5 <- eps[,5]
                   x.all <- data.frame(x1, x2, x3, x4, x5)
                   
                   out <- list()
                   for (s in 1:ncol(all.comb)){
                     # loop through settings
                     x <- x.all[,all.comb[,s]]
                     dat <- data.frame(y, x)
                     # assess well-specification
                     ms <- multi.spec(dat, B = n.split, return.predictor = FALSE, return.residual = FALSE,
                                      fitting = fitxg, predicting = predxg)
                     names(ms) <- paste(names(ms), paste(all.comb[,s], collapse=""), sep = "")
                     out <- c(out, ms)
                   }
                   out
                 }
    toc()
    stopCluster(cl)
    
    # initialize save file
    simulation <- list(n = n, r.seed = attr(res, "rng"), "commit" = commit)
    for (s in 1:ncol(all.comb)){
      # loop through settins
      mes <- paste(all.comb[,s], collapse="")
      # store output list to matrix
      res.steps.out <- array(unlist(res[, paste("steps", mes, sep = "")]),
                             dim = c(2 * n.split, p.mes, nsim), dimnames = list(NULL, colnames(res[1, paste("steps", mes, sep = "")][[1]]), NULL))
      res.sel.out <- array(unlist(res[, paste("sel", mes, sep = "")]), dim = c(2 * n.split, p.mes, nsim))
      res.pval <- matrix(unlist(res[, paste("pval", mes, sep = "")]), byrow = TRUE, nrow = nsim)
      res.pval.corr <- unlist(res[, paste("pval.corr", mes, sep = "")])
      simulation[[mes]] <- list(mes = all.comb[,s], steps.out = res.steps.out, sel.out = res.sel.out,
                             pval = res.pval, pval.corr = res.pval.corr)
      
      # print some summary results
      print(mes)
      print(median(res.pval.corr))
      for (j in 1:p.mes){
        print(paste(all.comb[j,s], ": ", sum(res.sel.out == j, na.rm = TRUE), sep = ""))
      }
    }
    
    
    # create unique filename based on sample size and time
    resname <- paste0("results n=", n, " ", format(Sys.time(), "%d-%b-%Y %H.%M"))
    # save the file to the folder
    if (save) save(simulation, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
  }
  return(paste("results/", newdir, sep = ""))
}