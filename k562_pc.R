library(reticulate)
require(parallel)
require(doRNG)
require(doSNOW)
require(git2r)
require(FOCI)
require(mgcv)
require(sfsmisc)
require(dHSIC)
require(xgboost)
source("multi_spec.R")
source("fitpred.R")
source("split.R")
np <- import("numpy")

table5 <- function(file, foc.mat = NULL){
  mat <- np$load(file) # load the preprocessed data
  obs <- which(mat[["interventions"]] == "non-targeting") # find observational data
  frac <- apply(mat[["expression_matrix"]][obs,] > 0, 2, mean) # active fraction per variables
  names <- mat[["var_names"]] # variable names
  all.env <- unique(mat["interventions"]) # all possible targets
  sel.col <- which(frac == 1) # restrict to always active variables
  n.col <- length(sel.col)
  sel.var <- names[sel.col]
  
  dat <- mat[["expression_matrix"]][obs, sel.col] # the data we work with
  colnames(dat) <- sel.var
  
  foci.var <- function(j){
    # estimate Markov blanket for the covariate
    focij <- foci(dat[,j], dat[,-j])
    sel <- focij$selectedVar$index
    sel[sel >= j] <- sel[sel >= j ] + 1
    # fix one of two options to be reproducible (seed does not help)
    if((j == 15 && (1 %in% sel)) || 
       (j == 24 & !(6 %in% sel)) ||
       (j == 7 & (2 %in% sel))){
      warning("Call again for reproducibility at j=", j)
      return(foci.var(j))
    }
    if (any(focij$selectedVar$names != sel.var[sel])) 
      stop("Non matching names")
    ret <- rep(0, ncol(dat))
    # encode in 0 one matrix
    ret[sel] <- 1
    return(ret)
  }
  
  if(is.null(foc.mat)){
    cat("Searching Markov boundaries")
    foc.mat <- foreach(j = 1:n.col, .combine = rbind) %do%{
      cat(paste(j, ""))
      foci.var(j)
    }
  } else if(!identical(dim(foc.mat), c(n.col, n.col))){
    stop("Bad matrix provided")
  }
  
  
  ad <- foc.mat + t(foc.mat) # check both directions, if it is in Markov boundary, this should be 2
  
  test.var <- function(i){
    # see if this variable affects others
    rows <- which(mat["interventions"] == sel.var[i]) # environment where this variable is intervened on
    dati <- mat[["expression_matrix"]][rows, sel.col]
    testi <- numeric(n.col)
    for (j in 1:n.col){
      testi[j] <- wilcox.test(dat[,j], dati[,j])$p.value # compare distributions
    }
    testi
  }
  cat("\n Mann-Whitney U tests \n")
  test <- foreach(j = 1:n.col, .combine = rbind) %do%{
    cat(paste(j, ""))
    test.var(j)
  }
  
  
  analyse.var <- function(y){
    subs <- sort(c(which(ad[y,] > 1), y)) # covariates which seem to be in Markov boundary
    test.sub <- test[subs, subs]
    
    if(length(subs) > 2 && any(test.sub[subs == y, subs != y] > 0.01)){ # not all are descendants
      diag(test.sub) <- 1
      
      dat.reg <- dat[,subs] # the data we work with
      wi <- which(colnames(dat.reg) == sel.var[y]) # find target position
      dat.reg <- cbind(dat.reg[,wi], dat.reg[,-wi]) # move target to first position
      colnames(dat.reg) <- c("y", colnames(dat[,subs])[-wi]) # call it y
      cat("Analysing target ", y, " with predictors ", paste(subs[-wi]), "\n")
      
      set.seed(4)
      # assess well-specification
      msx <- multi.spec(data.frame(dat.reg), B = 25, omit.global = FALSE, return.predictor = TRUE, return.residual = TRUE,
                        fitting = function(data, ind, res.ind) fitxg(data, ind, res.ind, ncol(data) - 1), predicting = predxg,
                        parallel = TRUE, sockets = 5, verbose = FALSE, export.functions = "fitxg")
      nsplit <- apply(!is.na(msx$steps), 2, sum) # how often is it not selected by FOCI
      pval.split <- fisher.split(is.na(msx$steps)) # proportion test
      pval <- test.sub[subs == y, subs != y] # Mann-Whitney-U test
      
      data.x <- xgb.DMatrix(dat.reg[,-1], label = dat.reg[,1])
      # estimate conditional mean on observational data
      fi.all <- xgb.cv(list(max_depth = ncol(dat.reg) - 1, nthread = 1), data = data.x, nrounds = 500, nfold = 2,
                       early_stopping_rounds = 3, prediction = TRUE, verbose = FALSE)
      
      env.data <- list()
      for (su in subs){
        # find the data where the predictors are intervened on
        rows <- which(mat["interventions"] == sel.var[su])
        dati <- mat[["expression_matrix"]][rows, sel.col[subs]]
        env.data[[sel.var[su]]] <- dati
      }
      
      watchlist <- env.data
      for (wa in names(watchlist)){
        # turn data to xgb compatible format
        watchlist[[wa]] <- xgb.DMatrix(watchlist[[wa]][,-which(subs == y)], label = watchlist[[wa]][,which(subs == y)])
      }
      # fit again for pre-defined rounds, evaluate on other environments
      fi.out <- xgb.train(list(max_depth = ncol(dat.reg) - 1, nthread = 1), data = data.x, nrounds = fi.all$best_iteration,
                          verbose = FALSE, watchlist = watchlist)
      
      par(mfrow = c(1 + 1 * (ncol(dat.reg) > 3), 2))
      bias <- c()
      for(wa in names(env.data)){
        if (wa == sel.var[y])
          next()
        pred <- predict(fi.out, watchlist[[wa]]) # predict on environment
        bias <- c(bias, mean((pred - env.data[[wa]][,which(subs == y)]))) # calculate bias
        plot(dat[,which(sel.var == wa)], dat[,y], xlim = c(0, 6), ylim = c(0, 6),
             xlab = wa, ylab = sel.var[y]) # observational data
        points(env.data[[wa]][,which(sel.var[subs] == wa)],  env.data[[wa]][,which(subs == y)], col = 2) # interventional data
        
        points(env.data[[wa]][,which(sel.var[subs] == wa)], pred, col = 3) # prediction
      }
      legend('bottomright', legend = c("observational", "interventional", "prediction"), col = 1:3, pch = 1)
      
      
      out <- cbind(pval, nsplit, pval.split, abs(bias)/mean(dat.reg[, 1])) # store to vector
      rownames(out) <- paste("$X_{", subs[subs != y], "}$", sep = "")
      # print LaTeX format
      cat("\\hline", "\n", sep = "")
      cat("\\multirow{", nrow(out), "}{*}{$X_{", y, "}$}", "\n", sep ="")
      for(i in 1:nrow(out)){
        cat("& ", rownames(out)[i], " & ", 
            paste(c(format(out[i,1], digits = 2, scientific = TRUE), out[i,2],
                    format(out[i,3], digits = 2, scientific = TRUE), format(out[i,4], digits = 2, scientific = TRUE)), collapse = " & "),
            " \\", "\\", "\n", sep = "")
      }
      out
    }  else {
      NULL
    }
  }
  
  cat("\n")
  # consider potential candidate y
  all.out <- list()
  for(y in 1:n.col){
    all.out[[y]] <- analyse.var(y)
  }
  
  y <- 0
  for (out in all.out){
    y <- y + 1
    if (is.null (out)) next()
    # print compactly in LaTeX format
    cat("\\hline", "\n", sep = "")
    cat("\\multirow{", nrow(out), "}{*}{$X_{", y, "}$}", "\n", sep ="")
    for(i in 1:nrow(out)){
      cat("& ", rownames(out)[i], " & ", 
          paste(c(format(out[i,1], digits = 2, scientific = TRUE), out[i,2],
                  format(out[i,3], digits = 2, scientific = TRUE), format(out[i,4], digits = 2, scientific = TRUE)), collapse = " & "),
          " \\", "\\", "\n", sep = "")
    }
  }
  return(list(all.out, foc.mat))
}


