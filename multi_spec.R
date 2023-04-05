multi.spec <- function(data, response = "y", B = 25, gamma = NULL, gamma.min = 0.05,
                       fitting = function(data) gam(wrapFormula(y ~., data = data), data = data)){
  if(is.null(gamma)){
    gamma <- (1 : (2 * B))/(2 * B)
    gamma <- gamma[gamma >= gamma.min]
  }
  n <- nrow(data)
  res.ind <- which(colnames(data) == response)
  pval <- numeric(2 * B)
  pred <- matrix(0, n, B)
  for (i in 1:B){
    cat((paste(i, " ")))
    ind <- sample(n, ceiling(n /2))
    fi1 <- fitting(data[ind,])
    fi2 <- fitting(data[-ind,])
    pred12 <- predict(fi1, newdata = data[-ind,])
    pred21 <- predict(fi2, newdata = data[ind,])
    pred[ind, i] <- pred21
    pred[-ind, i] <- pred12
    hs12 <- dhsic.test(data[-ind, res.ind] - pred12, data[-ind, -res.ind], method = "gamma")
    hs21 <- dhsic.test(data[ind, res.ind] - pred21, data[ind, -res.ind], method = "gamma")
    pval[c(i, i + B)] <- c(hs12$p.value, hs21$p.value)
  }
  pval
  quant.gamma <- quantile(pval, gamma, type = 1)/gamma
  penalty <- if (length(gamma) > 1) 
    (1 - log(min(gamma)))
  else 1
  pval.pre <- pmin(min(quant.gamma) * penalty, 1)
  # which.gamma <- which.min(quant.gamma)
  return(list(pval.corr = pval.pre, pval = pval, prediction = pred))
}