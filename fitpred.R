fitxg <- function(data, ind){
  n <- nrow(data)
  ind2 <- (1:n)[-ind]
  data.x <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1])
  xgb.cv(list(max_depth = 2, nthread = 1), data = data.x, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
}
predxg <- function(fit, data, ind){
  fit$pred[ind]
}

allfitxg <- function(data){
  data.x <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1])
  fi.all <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  out <- list()
  out$fitted.values <- fi.all$pred
  out$residuals <- data$y - fi.all$pred
  out
}

fitxg.het <- function(data, ind){
  n <- nrow(data)
  ind2 <- (1:n)[-ind]
  data.x <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1])
  fi1 <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  data.x2 <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1]^2)
  fi2 <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x2, nrounds = 500, folds = list(ind, ind2), early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  list(fi1, fi2)
}

pred.xg.het <- function(fit, data, ind){
  fit[[1]]$pred[ind]
}

norm.xg.het <- function(fit, data, ind){
  sqrt(fit[[2]]$pred[ind] - fit[[1]]$pred[ind]^2)
}

allfitxg.het <- function(data){
  data.x <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1])
  fi.all <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  data.x2 <- xgb.DMatrix(as.matrix(data[,-1]), label = data[,1]^2)
  fi.all2 <- xgb.cv(list(max_depth = 2, nthread = 1), data = data.x2, nrounds = 500, nfold = 2, early_stopping_rounds = 3, prediction = TRUE, verbose = F)
  out <- list()
  out$fitted.values <- fi.all$pred
  out$residuals <- (data$y - fi.all$pred)/sqrt(fi.all2$pred - fi.all$pred^2)
  out
}