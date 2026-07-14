# experiment #2
# question of interest: do different genes predict working memory differently across cortical brain network? 
# each cortical region is a node and edges represent structural connectivity
# response is average working memory activation

source("R/lasso.R")
source("R/rnc.R")
source("R/graph_basis.R")
source("R/sgvc.R")
source("R/network.lasso.R")
source("R/metrics.R")


X_df = read.csv("hcp-ahba/X_ahba_genes_complete_cases.csv",row.names=1,check.names=F)
A_df = read.csv("hcp-ahba/A_hcp_structural_ahba.csv",row.names=1,check.names=F)
y_df = read.csv("hcp-ahba/y_hcp_wm_2bk_minus_0bk.csv",row.names=1,check.names=F)

X = as.matrix(X_df)
A = as.matrix(A_df)
y = y_df$y_wm_2bk_0bk
#y = y_df$y_wm_2bk_0bk_z


n = nrow(X)
dinv = 1/sqrt(rowSums(A))
L = Diagonal(n)-Diagonal(x=dinv) %*%A%*% Diagonal(x=dinv)

set.seed(1)
test = sample(n, n%/%5)
train = setdiff(1:n, test)

t_sgvc.sparse = system.time(fit_sgvc.sparse <- sgvc(X, y, L, train, nb=20))[["elapsed"]]
t_lasso = system.time(fit_lasso <- lasso(X, y, train))[["elapsed"]]
t_enet = system.time(fit_enet <- elastic.net(X, y, train))[["elapsed"]]

### runtime (seconds)
rbind(sgvc.sparse=t_sgvc.sparse,lasso=t_lasso,enet=t_enet)


### out of sample results
rbind(sgvc.sparse=metrics(y[test],fit_sgvc.sparse$fitted[test]),
      lasso=metrics(y[test],fit_lasso$fitted[test]),
      enet=metrics(y[test],fit_enet$fitted[test])
)


fit_sgvc.sparse$nsel
fit_lasso$nsel
fit_enet$nsel

#fit_sgvc.sparse$sel_ind

