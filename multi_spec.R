multi.spec <- function(data, response = "y", B = 25, gamma = NULL, gamma.min = 0.05,
                       fitting = function(data, ind) gam(wrapFormula(y ~., data = data), data = data[ind, ]),
                       predicting = function(fit, data, ind) predict(fit, newdata = data[ind,]),
                       norming = NULL, trafo = abs, return.predictor = FALSE, return.residual = FALSE){
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
  for (i in 1:B){
    cat((paste(i, " ")))
    ind <- sample(n, ceiling(n /2))
    fi1 <- fitting(data, ind)
    fi2 <- fitting(data, (1:n)[-ind])
    pred12 <- predicting(fi1, data, -ind)
    pred21 <- predicting(fi2, data, ind)
    pred[ind, i] <- pred21
    pred[-ind, i] <- pred12
    
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
    
    res[ind, i] <- res21
    res[-ind, i] <- res12
    
    if (n /2 > 1e4){
      hs12 <- dhsic.test((res12)[1:1e4], data[-ind, -res.ind][1:1e4 ,], method = "gamma")
      hs21 <- dhsic.test((res21)[1:1e4], data[ind, -res.ind][1:1e4 ,], method = "gamma")
    } else {
      hs12 <- dhsic.test(res12, data[-ind, -res.ind], method = "gamma")
      hs21 <- dhsic.test(res21, data[ind, -res.ind], method = "gamma")
    }
    
    pval[c(i, i + B)] <- c(hs12$p.value, hs21$p.value)
    
    fo12 <- foci(trafo(res12), data[-ind, -res.ind])
    fo21 <- foci(trafo(res21), data[ind, -res.ind])
    steps12[i, fo12$selectedVar$index] <- diff(c(0, fo12$stepT))
    steps21[i, fo21$selectedVar$index] <- diff(c(0, fo21$stepT))
    sel12[i,] <- c(fo12$selectedVar$index, rep(NA, dim(data)[2] - 1 - length(fo12$selectedVar$index)))
    sel21[i,] <- c(fo21$selectedVar$index, rep(NA, dim(data)[2] - 1 - length(fo21$selectedVar$index)))
  }
  pval
  quant.gamma <- quantile(pval, gamma, type = 1)/gamma
  penalty <- if (length(gamma) > 1) 
    (1 - log(min(gamma)))
  else 1
  pval.pre <- pmin(min(quant.gamma) * penalty, 1)
  # which.gamma <- which.min(quant.gamma)
  out <- list(pval.corr = pval.pre, pval = pval,
              steps = rbind(steps21, steps12), sel = rbind(sel21, sel12))
  if (return.predictor) out$prediction <- pred
  if (return.residual) out$residual <- res
  return(out)
}