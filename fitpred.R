fitxg <- function(data, ind, res.ind, max_depth = 2){
  # function to fit first moment
  # Input
  # data (numeric, data frame): observational data
  # ind (integer, vector): indices to fit data on
  # res.ind (integer): column of target
  # max.depth (integer): tree depth
  # Output
  # xgb-fit
  n <- nrow(data)
  ind2 <- (1:n)[-ind]
  # transform data to xgb format
  data.x <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind])
  # fit data using predefined folds
  xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
}

predxg <- function(fit, data, ind){
  # function to extract prediction from xgb-fit
  fit$pred[ind]
}

allfitxg <- function(data, res.ind = 1, max_depth = 2, verbose = FALSE){
  # function to fit first moment on all data
  # Input
  # data (numeric, data frame): observational data
  # res.ind (integer): column of target
  # max.depth (integer): tree depth
  # verbose (boolean): print xgb information
  # Output
  # fitted.values (numeric, vector): estimated first moment
  # residuals (numeric, vector): estimated residuals
  data.x <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind])
  fi.all <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = verbose)
  out <- list()
  out$fitted.values <- fi.all$pred
  out$residuals <- data[, res.ind] - fi.all$pred
  out
}

fitxg.het <- function(data, ind, res.ind, max_depth = 2){
  # function to fit first and second moment
  # Input
  # data (numeric, data frame): observational data
  # ind (integer, vector): indices to fit data on
  # res.ind (integer): column of target
  # max.depth (integer): tree depth
  # Output
  # xgb-fits for both moments
  n <- nrow(data)
  ind2 <- (1:n)[-ind]
  data.x <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind])
  fi1 <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  data.x2 <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind]^2)
  fi2 <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x2, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  list(fi1, fi2)
}

predxg.het <- function(fit, data, ind){
  # function to extract prediction from xgb-fit
  fit[[1]]$pred[ind]
}

normxg.het <- function(fit, data, ind){
  # function to extract estimated standard deviation from xgb-fits
  sqrt(fit[[2]]$pred[ind] - fit[[1]]$pred[ind]^2)
}

allfitxg.het <- function(data, res.ind = 1, max_depth = 2, verbose = FALSE){
  # function to fit first and second moment on all data and get normalized residuals
  # Input
  # data (numeric, data frame): observational data
  # res.ind (integer): column of target
  # max.depth (integer): tree depth
  # verbose (boolean): print xgb information
  # Output
  # fitted.values (numeric, vector): estimated first moment
  # residuals (numeric, vector): estimated residuals
  data.x <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind])
  fi.all <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = verbose)
  data.x2 <- xgb.DMatrix(as.matrix(data[,-res.ind]), label = data[,res.ind]^2)
  fi.all2 <- xgb.cv(list(max_depth = max_depth, nthread = 1), data = data.x2, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = verbose)
  out <- list()
  out$fitted.values <- fi.all$pred
  resu <- data[, res.ind] - fi.all$pred
  norm <- sqrt(fi.all2$pred - fi.all$pred^2)
  out$residuals <- resu/norm
  if (any(is.na(out$residuals)))
    out$residuals[is.na(out$residuals)] <- 10 * max(abs(out$residuals), na.rm = TRUE) * sign(resu[is.na(out$residuals)])
  out
}