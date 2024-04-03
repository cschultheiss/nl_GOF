multi.spec <- function(data, response = "y", B = 25, gamma = NULL, gamma.min = 0.05,
                       fitting = function(data, ind, res.ind = NULL) gam(wrapFormula(y ~., data = data), data = data[ind, ]),
                       predicting = function(fit, data, ind) predict(fit, newdata = data[ind,]),
                       norming = NULL, omit.global = FALSE, trafo = abs, return.predictor = FALSE, 
                       return.residual = FALSE, return.indices = FALSE, verbose = TRUE, parallel = FALSE, sockets = NULL, 
                       export.packages = c("mgcv", "sfsmisc", "xgboost", "FOCI", "dHSIC")){
  
  # function to assess well-specification
  # Input
  # data (numeric, data frame): observational data, requires named columns
  # response (string): column name of response variable
  # B (integer): number of splits
  # gamma (numeric, vector): quantiles to be considered for global test
  # gamma.min (numeric): minimum quantile to be considered for global test
  # fitting (function): estimator of conditional mean and potentially second moment
  # predicting (function): predictor on other split based on fit
  # norming (function): estimator of conditional standard deviation
  # omit.global (boolean): omit global test?
  # trafo (function): transformation on (estimated) residuals
  # return.predictor (boolean): return the estimated conditional mean?
  # return.residual (boolean): return the estimated residuals?
  # return.indices (boolean): return the random splits?
  # verbose (boolean): print intermediate information?
  # parallel (boolean): apply parralelization?
  # sockets (integer): number of sockets for parallelization
  # export.packages (character, vector): packages needed in sockets
  # Output
  # pval.corr (numeric): combined p-value
  # pval (numeric, vector): individual p-values
  # steps (numeric, matrix): stepsize corresponding to that variable
  # sel (integer, matrix): order of variable selection
  
  if(parallel && is.null(sockets))
    stop("Need to provide number of sockets for parallelization")
  if(is.null(gamma)){
    # use all useful quantiles if none defined
    gamma <- (1 : (2 * B))/(2 * B)
    gamma <- gamma[gamma >= gamma.min]
  }
  
  n <- nrow(data) # number of observations
  res.ind <- which(colnames(data) == response) # where to find response
  pval <- numeric(2 * B)
  res <- pred <- matrix(0, n, B) # residuals and predictors
  sel12 <- sel21 <- matrix(NA, B, dim(data)[2] - 1) # store if variable is selected
  steps12 <-  matrix(NA, nrow = B, ncol = dim(data)[2] - 1) # store step size
  colnames(steps12) <- colnames(data)[-res.ind]
  steps21 <- steps12
  
  # function for one split to be repeated
  one.split <- function(){
    ind <- sample(n, ceiling(n /2)) # random split
    fi1 <- fitting(data, ind, res.ind) # fit on first half
    fi2 <- fitting(data, (1:n)[-ind], res.ind) # fit on second half
    pred12 <- predicting(fi1, data, -ind) # predict on the other
    pred21 <- predicting(fi2, data, ind)
    res12 <- data[-ind, res.ind] - pred12 # get residuum
    res21 <- data[ind, res.ind] - pred21
    
    if(!is.null(norming)){
      norm12 <- norming(fi1, data, -ind) # get conditional standard deviation
      norm21 <- norming(fi2, data, ind)
      res12u <- res12
      res21u <- res21
      res12 <- res12 / norm12 # normalize residual
      # replace NAs by large value
      if (any(is.na(res12)))
        res12[is.na(res12)] <- 10 * max(abs(res12), na.rm = TRUE) * sign(res12u[is.na(res12)])
      res21  <- res21 / norm21
      if (any(is.na(res21)))
        res21[is.na(res21)] <- 10 * max(abs(res21), na.rm = TRUE) * sign(res21u[is.na(res21)])
    }
    
    if(!omit.global){
      if (n /2 > 1e4){ # use subdata if it is too large
        # global test
        hs12 <- dhsic.test((res12)[1:1e4], data[-ind, -res.ind][1:1e4 ,], method = "gamma")
        hs21 <- dhsic.test((res21)[1:1e4], data[ind, -res.ind][1:1e4 ,], method = "gamma")
      } else {
        hs12 <- dhsic.test(res12, data[-ind, -res.ind], method = "gamma")
        hs21 <- dhsic.test(res21, data[ind, -res.ind], method = "gamma")
      }
    } else {
      hs12 <- hs21 <- list(p.value = 1)
    }
    
    # apply FOCI
    fo12 <- foci(trafo(res12), subset(data[-ind,], select = - res.ind))
    fo21 <- foci(trafo(res21), subset(data[ind,], select = - res.ind))
    
    # get step size
    steps12 <- steps21 <-  rep(NA, dim(data)[2] - 1)
    steps12[fo12$selectedVar$index] <- diff(c(0, fo12$stepT))
    steps21[fo21$selectedVar$index] <- diff(c(0, fo21$stepT))
    
    out <- list(pval12 = hs12$p.value, pval21 = hs21$p.value,
                steps12 = steps12, steps21 = steps21,
                sel12 = c(fo12$selectedVar$index, rep(NA, dim(data)[2] - 1 - length(fo12$selectedVar$index))),
                sel21 = c(fo21$selectedVar$index, rep(NA, dim(data)[2] - 1 - length(fo21$selectedVar$index))))
    if (return.predictor){
      pred <- numeric(n)
      pred[ind] <- pred21
      pred[-ind] <- pred12
      out$pred <- pred
    }
    if (return.residual){
      res <- numeric(n)
      res[ind] <- res21
      res[-ind] <- res12
      out$res <- res
    }
    if(return.indices)
      out$ind <- ind
    out
  }
  
  # apply multiple times
  if (parallel) {
    progress <- function(n, tag) {
      if(verbose) cat(paste(n , ""))
    }
    
    opts <- list(progress = progress)
    cl<-makeSOCKcluster(sockets)
    registerDoSNOW(cl)
    out.splits <- foreach(i = 1:B, .combine = rbind,
                          .packages = export.packages, .options.snow = opts) %dorng%{
      one.split()
    }
    stopCluster(cl)
  } else {
    out.splits <- foreach(i = 1:B, .combine = rbind) %do%{
      if(verbose) cat(paste(i, ""))
      one.split()
    }
  }
    
  # combined p-value
  pval <- c(unlist(out.splits[,"pval12"]), unlist(out.splits[,"pval21"]), use.names = FALSE) 
  quant.gamma <- quantile(pval, gamma, type = 1)/gamma
  penalty <- if (length(gamma) > 1) 
    (1 - log(min(gamma)))
  else 1
  pval.pre <- pmin(min(quant.gamma) * penalty, 1)
  
  # store to matrices
  steps12 <- matrix(unlist(out.splits[, "steps12"]), byrow = TRUE, nrow = B, dimnames = list(NULL, colnames(data)[-res.ind]))
  steps21 <- matrix(unlist(out.splits[, "steps21"]), byrow = TRUE, nrow = B, dimnames = list(NULL, colnames(data)[-res.ind]))
  sel12 <- matrix(unlist(out.splits[, "sel12"]), byrow = TRUE, nrow = B)
  sel21 <- matrix(unlist(out.splits[, "sel21"]), byrow = TRUE, nrow = B)
  
  
  out <- list(pval.corr = pval.pre, pval = pval,
              steps = rbind(steps21, steps12), sel = rbind(sel21, sel12))
  if (return.predictor) out$prediction <- matrix(unlist(out.splits[, "pred"]), byrow = FALSE, nrow = n)
  if (return.residual) out$residual <- matrix(unlist(out.splits[, "res"]), byrow = FALSE, nrow = n)
  if (return.indices) out$indices <- matrix(unlist(out.splits[, "ind"]), byrow = FALSE, nrow = ceiling(n /2))
  return(out)
}