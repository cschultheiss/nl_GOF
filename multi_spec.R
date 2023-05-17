multi.spec <- function(data, response = "y", B = 25, gamma = NULL, gamma.min = 0.05,
                       fitting = function(data, data.val) gam(wrapFormula(y ~., data = data), data = data),
                       predicting = function(fit, data.val) predict(fit, newdata = data.val),
                       return.predictor = FALSE){
  if(is.null(gamma)){
    gamma <- (1 : (2 * B))/(2 * B)
    gamma <- gamma[gamma >= gamma.min]
  }
  n <- nrow(data)
  res.ind <- which(colnames(data) == response)
  pval <- numeric(2 * B)
  pred <- matrix(0, n, B)
  sel12 <- sel21 <- matrix(NA, B, dim(data)[2] - 1)
  steps12 <-  matrix(NA, nrow = B, ncol = dim(data)[2] - 1)
  colnames(steps12) <- colnames(data)[-res.ind]
  steps21 <- steps12
  for (i in 1:B){
    cat((paste(i, " ")))
    ind <- sample(n, ceiling(n /2))
    fi1 <- fitting(data[ind,], dat[-ind,])
    fi2 <- fitting(data[-ind,], data[ind,])
    pred12 <- predicting(fi1, data.val = data[-ind,])
    pred21 <- predicting(fi2, data.val = data[ind,])
    pred[ind, i] <- pred21
    pred[-ind, i] <- pred12
    
    if (n /2 > 1e4){
      hs12 <- dhsic.test((data[-ind, res.ind] - pred12)[1:1e4], data[-ind, -res.ind][1:1e4 ,], method = "gamma")
      hs21 <- dhsic.test((data[ind, res.ind] - pred21)[1:1e4], data[ind, -res.ind][1:1e4 ,], method = "gamma")
    } else {
      hs12 <- dhsic.test(data[-ind, res.ind] - pred12, data[-ind, -res.ind], method = "gamma")
      hs21 <- dhsic.test(data[ind, res.ind] - pred21, data[ind, -res.ind], method = "gamma")
    }
    
    pval[c(i, i + B)] <- c(hs12$p.value, hs21$p.value)
    
    fo12 <- foci(abs(data[-ind, res.ind] - pred12), data[-ind, -res.ind])
    fo21 <- foci(abs(data[ind, res.ind] - pred21), data[ind, -res.ind])
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
  return(out)
}