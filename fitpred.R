fitxg <- function(data, ind, res.ind, max_depth = 2){
  n <- nrow(data)
  ind2 <- (1:n)[-ind]
  data.x <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind])
  xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
}
predxg <- function(fit, data, ind){
  fit$pred[ind]
}

allfitxg <- function(data, max_depth = 2, verbose = FALSE){
  data.x <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1])
  fi.all <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = verbose)
  out <- list()
  out$fitted.values <- fi.all$pred
  out$residuals <- data$y - fi.all$pred
  out
}

fitxg.het <- function(data, ind, res.ind, max_depth = 2){
  n <- nrow(data)
  ind2 <- (1:n)[-ind]
  data.x <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind])
  fi1 <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  data.x2 <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind]^2)
  fi2 <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x2, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  list(fi1, fi2)
}

predxg.het <- function(fit, data, ind){
  fit[[1]]$pred[ind]
}

normxg.het <- function(fit, data, ind){
  sqrt(fit[[2]]$pred[ind] - fit[[1]]$pred[ind]^2)
}

allfitxg.het <- function(data){
  data.x <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1])
  fi.all <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  data.x2 <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1]^2)
  fi.all2 <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x2, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  out <- list()
  out$fitted.values <- fi.all$pred
  resu <- data$y - fi.all$pred
  norm <- sqrt(fi.all2$pred - fi.all$pred^2)
  out$residuals <- resu/norm
  if (any(is.na(out$residuals)))
    out$residuals[is.na(out$residuals)] <- 10 * max(abs(out$residuals), na.rm = TRUE) * sign(resu[is.na(out$residuals)])
  out
}