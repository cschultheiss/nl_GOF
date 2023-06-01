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