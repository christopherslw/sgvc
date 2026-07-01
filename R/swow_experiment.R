# experiment #1
# question of interest: do lexical features predict reaction time differently across the SWOW semantic network? 
# each word is a node and edges represent strong human free associations 
# response is lexical decision reaction time

source("lasso.R")
source("rnc.R")
source("graph_basis.R")
source("sgvc.sparse.R")
source("sgvc.ls.R")
source("sgvc.ls.cv.R")
source("network.lasso.R")
source("cv_out.R")
source("cv_nb.R")
source("metrics.R")


node = read.csv("swow/swow_elp_reaction_time_node_data.csv")
node$z_lgsubtlwf[is.na(node$z_lgsubtlwf)] = 0   # only one missing value -> impute
feats = c("z_length", "z_log_freq_hal", "z_lgsubtlwf", "z_ortho_n", "z_old", "z_nphon", "z_nsyll", "z_nmorph")
X = as.matrix(node[, feats])
y = node$y_ldt_z
n = nrow(X)
edge = read.csv("swow/swow_elp_reaction_time_edge_list.csv")
A = matrix(0, n, n)
A[cbind(edge$from, edge$to)] = edge$weight
A[cbind(edge$to, edge$from)] = edge$weight


nb = 9
nb_max = min(20, n-2) # search range for the spectral basis length
dinv = 1 / sqrt(rowSums(A))
L = diag(n) - A*outer(dinv, dinv)
Phi = graph_basis(L, nb_max)
set.seed(1)
fold0 = sample(rep(1:5, length.out=n)) # choose the spectral basis length nb by cv
t_ls = system.time(nb_ls <- sgvc.ls.cv(X, y, L, fold0, 1:nb_max))[["elapsed"]]
t_sparse = system.time(nb_sparse <- cv_nb(function(tr, nb) sgvc.sparse(X, y, L, tr, nb), 1:nb_max, fold0, y))[["elapsed"]]


# cross validation on different methods

methods = list(lasso=function(tr) lasso(X, y, tr),
               rnc=function(tr) rnc(X, y, Phi[, 1:nb], tr),
               sgvc.ls=function(tr) sgvc.ls(X, y, L, tr, nb_ls),
               sgvc.sparse=function(tr) sgvc.sparse(X, y, L, tr, nb_sparse),
               network.lasso = function(tr) network.lasso(X, y, A, tr))
sel = c(lasso=0, rnc=0, sgvc.ls=t_ls, sgvc.sparse=t_sparse, network.lasso=0) # factor in time for cv_nb

reps = 5
set.seed(2)
folds = lapply(1:reps, function(r) sample(rep(1:5, length.out=n)))
score = array(0, dim=c(length(methods), 3, reps))
time = numeric(length(methods))

for (mthd in 1:length(methods)) {
  t = system.time(
    for (r in 1:reps) {
      yhat = cv_out(methods[[mthd]], folds[[r]])
      score[mthd,,r] = metrics(yhat)
    }
  )[["elapsed"]]
  time[mthd] = t + sel[[names(methods)[mthd]]]
}

tab = apply(score, c(1, 2), mean)
rownames(tab) = names(methods)
colnames(tab) = c("R2", "RMSE", "MAE")
tab = cbind(tab, seconds=time)
print(round(tab[order(-tab[,"R2"]),], 3))
