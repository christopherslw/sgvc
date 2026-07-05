# experiment #1
# question of interest: do lexical features predict reaction time differently across the SWOW semantic network? 
# each word is a node and edges represent strong human free associations 
# response is lexical decision reaction time
source("R/lasso.R")
source("R/rnc.R")
source("R/graph_basis.R")
source("R/sgvc.R")
source("R/network.lasso.R")
source("R/metrics.R")


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

dinv = 1/sqrt(rowSums(A))
L = Diagonal(n)-Diagonal(x=dinv) %*%A%*% Diagonal(x=dinv)

set.seed(1)
test = sample(n, n%/%5)
train = setdiff(1:n, test)

t_sgvc.ls = system.time(fit_sgvc.ls <- sgvc.ls(X, y, L, train))[["elapsed"]]
t_sgvc.sparse = system.time(fit_sgvc.sparse <- sgvc(X, y, L, train, nb=fit_sgvc.ls$nb))[["elapsed"]]
t_lasso = system.time(fit_lasso <- lasso(X, y, train))[["elapsed"]]
t_enet = system.time(fit_enet <- elastic.net(X, y, train))[["elapsed"]]
t_netlasso = system.time(fit_netlasso <- network.lasso(X, y, A, train))[["elapsed"]]
t_rnc = system.time(fit_rnc <- rnc(X, y, L, train))[["elapsed"]]


### runtime (seconds)
rbind(sgvc.ls=t_sgvc.ls,sgvc.sparse=t_sgvc.sparse,lasso=t_lasso,enet=t_enet,netlasso=t_netlasso,rnc=t_rnc)


### out of sample results
rbind(sgvc.ls=metrics(y[test],fit_sgvc.ls$fitted[test]),
      sgvc.sparse=metrics(y[test],fit_sgvc.sparse$fitted[test]),
      lasso=metrics(y[test],fit_lasso$fitted[test]),
      enet=metrics(y[test],fit_enet$fitted[test]),
      netlasso=metrics(y[test],fit_netlasso$fitted[test]),
      rnc=metrics(y[test],fit_rnc$fitted[test])
)

