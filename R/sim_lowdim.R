library(Matrix)
source("R/lasso.R")
source("R/elastic.net.R")
source("R/graph_basis.R")
source("R/sgvc.R")
source("R/gen_data.R")
source("R/metrics.R")
source("R/network.lasso.R")
source("R/rnc.R")

set.seed(1) 
n = 500
p = 6
s = 6  # num of true active covariates
data = gen_data(n, p, s, graph="sbm")
X = data$X
y = data$y
n = data$n
A = data$W

# normalized Laplacian kept sparse
dinv = 1/sqrt(rowSums(A))
L = Diagonal(n)-Diagonal(x=dinv) %*%A%*% Diagonal(x=dinv)

test = sample(n, n%/%5)
train = setdiff(1:n, test)

#t_sgvc.ls = system.time(fit_sgvc.ls <- sgvc.ls(X, y, L, train))[["elapsed"]]
t_sgvc.sparse = system.time(fit_sgvc.sparse <- sgvc(X, y, L, train))[["elapsed"]]
t_ols = system.time(fit_ols <- ols(X, y, train))[["elapsed"]]
t_lasso = system.time(fit_lasso <- lasso(X, y, train))[["elapsed"]]
t_enet = system.time(fit_enet <- elastic.net(X, y, train))[["elapsed"]]
t_netlasso = system.time(fit_netlasso <- network.lasso(X, y, A, train))[["elapsed"]]
t_rnc = system.time(fit_rnc <- rnc(X, y, L, train))[["elapsed"]]


### in sample results

# R2,RMSE,MAE
rbind(sgvc.sparse=metrics(y[train],fit_sgvc.sparse$fitted[train]),
      ols=metrics(y[train],fit_ols$fitted[train]),
      lasso=metrics(y[train],fit_lasso$fitted[train]),
      enet=metrics(y[train],fit_enet$fitted[train]),
      netlasso=metrics(y[train],fit_netlasso$fitted[train]),
      rnc=metrics(y[train],fit_rnc$fitted[train])
)
#runtime (seconds)
rbind(sgvc.sparse=t_sgvc.sparse,ols=t_ols,lasso=t_lasso,enet=t_enet,netlasso=t_netlasso,rnc=t_rnc)


### out of sample results

# R2,RMSE,MAE
rbind(sgvc.sparse=metrics(y[test],fit_sgvc.sparse$fitted[test]),
      ols=metrics(y[test],fit_ols$fitted[test]),
      lasso=metrics(y[test],fit_lasso$fitted[test]),
      enet=metrics(y[test],fit_enet$fitted[test]),
      netlasso=metrics(y[test],fit_netlasso$fitted[test]),
      rnc=metrics(y[test],fit_rnc$fitted[test])
)


