multi.spec <- function(data, response = "y", B = 25, gamma = NULL, gamma.min = 0.05,
                       fitting = function(data, ind) gam(wrapFormula(y ~., data = data), data = data[ind, ]),
                       predicting = function(fit, data, ind) predict(fit, newdata = data[ind,]),
                       norming = NULL, trafo = abs, return.predictor = FALSE, return.residual = FALSE,
                       return.indices =FALSE, parallel = FALSE, sockets = NULL){
  if(parallel && is.null(sockets))
    stop("Need to provide number of sockets for parallelization")
  if(is.null(gamma)){
    gamma <- (1 : (2 * B))/(2 * B)
    gamma <- gamma[gamma >= gamma.min]
  }
  n <- nrow(data)
  res.ind <- which(colnames(data) == response)
  pval <- numeric(2 * B)
  res <- pred <- matrix(0, n, B)
  sel12 <- sel21 <- matrix(NA, B, dim(data)[2] - 1)
  steps12 <-  matrix(NA, nrow = B, ncol = dim(data)[2] - 1)
  colnames(steps12) <- colnames(data)[-res.ind]
  steps21 <- steps12
  
  one.split <- function(){
    ind <- sample(n, ceiling(n /2))
    fi1 <- fitting(data, ind)
    fi2 <- fitting(data, (1:n)[-ind])
    pred12 <- predicting(fi1, data, -ind)
    pred21 <- predicting(fi2, data, ind)
    res12 <- data[-ind, res.ind] - pred12
    res21 <- data[ind, res.ind] - pred21
    
    if(!is.null(norming)){
      norm12 <- norming(fi1, data, -ind)
      norm21 <- norming(fi2, data, ind)
      res12u <- res12
      res21u <- res21
      res12 <- res12 / norm12
      if (any(is.na(res12)))
        res12[is.na(res12)] <- 10 * max(abs(res12), na.rm = TRUE) * sign(res12u[is.na(res12)])
      res21  <- res21 / norm21
      if (any(is.na(res21)))
        res21[is.na(res21)] <- 10 * max(abs(res21), na.rm = TRUE) * sign(res21u[is.na(res21)])
    }
    if (n /2 > 1e4){
      hs12 <- dhsic.test((res12)[1:1e4], data[-ind, -res.ind][1:1e4 ,], method = "gamma")
      hs21 <- dhsic.test((res21)[1:1e4], data[ind, -res.ind][1:1e4 ,], method = "gamma")
    } else {
      hs12 <- dhsic.test(res12, data[-ind, -res.ind], method = "gamma")
      hs21 <- dhsic.test(res21, data[ind, -res.ind], method = "gamma")
    }
    
    fo12 <- foci(trafo(res12), subset(data[-ind,], select = - res.ind))
    fo21 <- foci(trafo(res21), subset(data[ind,], select = - res.ind))
    
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
  
  if (parallel) {
    progress <- function(n, tag) {
      cat(paste(n , ""))
    }
    
    opts <- list(progress = progress)
    cl<-makeSOCKcluster(sockets)
    registerDoSNOW(cl)
    out.splits <- foreach(i = 1:B, .combine = rbind,
                          .packages = c("mgcv", "sfsmisc", "xgboost", "FOCI", "dHSIC"), .options.snow = opts) %dorng%{
      one.split()
    }
    stopCluster(cl)
  } else {
    out.splits <- foreach(i = 1:B, .combine = rbind) %do%{
      cat(paste(i, ""))
      one.split()
    }
  }
    
  pval <- c(unlist(out.splits[,"pval12"]), unlist(out.splits[,"pval21"]), use.names = FALSE) 
  quant.gamma <- quantile(pval, gamma, type = 1)/gamma
  penalty <- if (length(gamma) > 1) 
    (1 - log(min(gamma)))
  else 1
  pval.pre <- pmin(min(quant.gamma) * penalty, 1)
  # which.gamma <- which.min(quant.gamma)
  
  steps12 <- matrix(unlist(out.splits[, "steps12"]), byrow = TRUE, nrow = B, dimnames = list(NULL, colnames(data)[-1]))
  steps21 <- matrix(unlist(out.splits[, "steps21"]), byrow = TRUE, nrow = B, dimnames = list(NULL, colnames(data)[-1]))
  sel12 <- matrix(unlist(out.splits[, "sel12"]), byrow = TRUE, nrow = B)
  sel21 <- matrix(unlist(out.splits[, "sel21"]), byrow = TRUE, nrow = B)
  
  
  out <- list(pval.corr = pval.pre, pval = pval,
              steps = rbind(steps21, steps12), sel = rbind(sel21, sel12))
  if (return.predictor) out$prediction <- matrix(unlist(out.splits[, "pred"]), byrow = FALSE, nrow = n)
  if (return.residual) out$residual <- matrix(unlist(out.splits[, "res"]), byrow = FALSE, nrow = n)
  if (return.indices) out$indices <- matrix(unlist(out.splits[, "ind"]), byrow = FALSE, nrow = ceiling(n /2))
  return(out)
}