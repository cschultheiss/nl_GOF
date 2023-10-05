library(reticulate)
require(parallel)
require(pcalg)
require(kpcalg)
require(ParallelPC)
require(git2r)
require(tictoc)
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

mat <- np$load("data/dataset_k562_filtered.npz")

obs <- which(mat[["interventions"]] == "non-targeting") 

frac <- apply(mat[["expression_matrix"]][obs,] > 0, 2, mean)

frac.all <- apply(mat[["expression_matrix"]] > 0, 2, mean)

names <- mat[["var_names"]]

all.env <- unique(mat["interventions"])

sel.col <- which(frac == 1)
n.col <- length(sel.col)
sel.var <- names[sel.col]

dat <- mat[["expression_matrix"]][obs, sel.col]
colnames(dat) <- sel.var


# tic()
# algo <- rfci_parallel(suffStat = list(data=dat, ic.method="hsic.gamma"),
#    indepTest = kernelCItest, alpha = 1e-3, m.max = 2,
#    labels = sel.var, verbose = TRUE, num.cores = 20)
# toc()

# save <- TRUE
# # create save location, adjust depending on folder structure
# if (save) {
#   commit <- revparse_single(revision = "HEAD")
#   newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
#   dir.create(paste("results/", newdir, sep=""))
#   resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"))
#   out <- list(algo = algo, commit = commit)
#   save(out, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
# }

foc.mat <- matrix(0, n.col, n.col)
for (j in 1:n.col){
  print(j)
  focij <- foci(dat[,j], dat[,-j])
  sel <- focij$selectedVar$index
  sel[sel >= j] <- sel[sel >= j ] + 1
  foc.mat[j, sel] <- 1
  if (any(focij$selectedVar$names != sel.var[sel])) 
    stop("Non matching names")
}

ad <- foc.mat + t(foc.mat)

test <- matrix(NA, n.col, n.col)
for(i in 1:n.col){
  print(i)
  rows <- which(mat["interventions"] == sel.var[i])
  dati <- mat[["expression_matrix"]][rows,sel.col]
  for (j in 1:n.col){
   test[i,j] <- wilcox.test(dat[,j], dati[,j])$p.value
  }
  print(log10(test[i,]))
}


# consider potential candidate y
for(y in 1:28){
  subs <- sort(c(which(ad[y,] > 1),y))
  
  test.sub <- test[subs, subs]
  if(length(subs) > 1)
    cat(y, ": ", test.sub[subs == y, subs != y], "\n")
}

y <- 16 # 6, 9 also work "reasonable"
subs <- sort(c(which(ad[y,] > 1),y))

test.sub <- test[subs, subs]

diag(test.sub) <- 1


dat.reg <- dat[,subs]
wi <- which(colnames(dat.reg) == sel.var[y])
dat.reg <- cbind(dat.reg[,wi], dat.reg[,-wi])
colnames(dat.reg) <- c("y", colnames(dat[,subs])[-wi])

set.seed(4)
msx <- multi.spec(data.frame(dat.reg), B = 25, return.predictor = TRUE, return.residual = TRUE,
  fitting = function(data, ind) fitxg(data, ind, ncol(data) - 1), predicting = predxg, parallel = TRUE, sockets = 5)
apply(is.na(msx$steps), 2, sum)
fisher.split(is.na(msx$steps))
mean(msx$residual^2)
var(c(msx$residual))

msg <- multi.spec(data.frame(dat.reg), B = 25, return.predictor = TRUE, return.residual = TRUE,
                  parallel = TRUE, sockets = 5)
apply(is.na(msg$steps), 2, sum)
fisher.split(is.na(msg$steps))
mean(msg$residual^2)
var(c(msg$residual))


data.x <- xgb.DMatrix(dat.reg[,-1], label = dat.reg[,1])
fi.all <- xgb.cv(list(max_depth = ncol(dat.reg) - 1, nthread = 1), data = data.x, nrounds = 500, nfold = 2,
                 early_stopping_rounds = 3, prediction = TRUE, verbose = TRUE)



env.data <- list()
for (su in subs){
  rows <- which(mat["interventions"] == sel.var[su])
  dati <- mat[["expression_matrix"]][rows, sel.col[subs]]
  env.data[[sel.var[su]]] <- dati
}

watchlist <- env.data
for (wa in names(watchlist)){
  watchlist[[wa]] <- xgb.DMatrix(watchlist[[wa]][,-which(subs == y)], label = watchlist[[wa]][,which(subs == y)])
}
fi.out <- xgb.train(list(max_depth = ncol(dat.reg) - 1, nthread = 1), data = data.x, nrounds = fi.all$best_iteration,
                    verbose = TRUE, watchlist = watchlist)

par(mfrow = c(3,1))
for(wa in names(env.data)){
  if (wa == sel.var[y])
    next()
  pred <- predict(fi.out, watchlist[[wa]])
  print(mean((pred - env.data[[wa]][,which(subs == y)])))

  # qqplot(dat[,y], env.data[[wa]][,which(subs == y)])
  # abline(0,1)
  # qqplot(dat[,which(sel.var == wa)], env.data[[sel.var[y]]][,which(sel.var[subs] == wa)])
  # abline(0,1)
  plot(dat[,which(sel.var == wa)], dat[,y], xlim = c(0, 6), ylim = c(0, 6),
       xlab = wa, ylab = sel.var[y])
  points(env.data[[wa]][,which(sel.var[subs] == wa)],  env.data[[wa]][,which(subs == y)], col = 2)
  # points(env.data[[sel.var[y]]][,which(sel.var[subs] == wa)],  env.data[[sel.var[y]]][,which(subs == y)], col = 3)
  points(env.data[[wa]][,which(sel.var[subs] == wa)], pred, col = 3)
}

