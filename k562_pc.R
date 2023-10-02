library(reticulate)
require(parallel)
require(pcalg)
require(kpcalg)
require(ParallelPC)
require(git2r)
require(tictoc)

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

tic()
algo <- rfci_parallel(suffStat = list(data=dat, ic.method="hsic.gamma"),
   indepTest = kernelCItest, alpha = 1e-3, m.max = 2,
   labels = sel.var, verbose = TRUE, num.cores = 20)
toc()

save <- TRUE
# create save location, adjust depending on folder structure
if (save) {
  commit <- revparse_single(revision = "HEAD")
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("results/", newdir, sep=""))
  resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"))
  out <- list(algo = algo, commit = commit)
  save(out, file = paste("results/", newdir, "/", resname, ".RData", sep = ""))
}

# test <- matrix(NA, n.col, n.col)
# for(i in 1:n.col){
#   print(i)
#   rows <- which(mat["interventions"] == sel.var[i])
#   dati <- mat[["expression_matrix"]][rows,sel.col]
#   for (j in 1:n.col){
#    test[i,j] <- wilcox.test(dat[,j], dati[,j])$p.value
#   }
#   print(log10(test[i,]))
# }
