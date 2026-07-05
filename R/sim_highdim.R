library(Matrix)
source("R/lasso.R")
source("R/elastic.net.R")
source("R/graph_basis.R")
source("R/sgvc.R")
source("R/gen_data.R")
source("R/metrics.R")


set.seed(1) 
n = 5000
p = 20000
s = 12  # num of true active covariates
data = gen_data(n, p, s, graph="lattice")
X = data$X
y = data$y
n = data$n
A = data$W

# normalized Laplacian kept sparse
dinv = 1/sqrt(rowSums(A))
L = Diagonal(n)-Diagonal(x=dinv) %*%A%*% Diagonal(x=dinv)

test = sample(n, n%/%5)
train = setdiff(1:n, test)

t_sgvc = system.time(fit_sgvc.sparse <- sgvc(X, y, L, train))[["elapsed"]]
t_lasso = system.time(fit_lasso <- lasso(X, y, train))[["elapsed"]]
t_enet = system.time(fit_enet <- elastic.net(X, y, train))[["elapsed"]]



### runtime (seconds)
rbind(sgvc.sparse=t_sgvc,lasso=t_lasso,enet=t_enet)


### out of sample prediction and variable selection 
sel_rates = function(sel, truth, p) {
  tpr = mean(truth %in% sel)
  fpr = length(setdiff(sel, truth)) / (p-length(truth))
  c(TPR=tpr, FPR=fpr, nsel=length(sel))
}
truth = data$active
#R2,RMSE,MAE
rbind(sgvc.sparse=metrics(y[test],fit_sgvc.sparse$fitted[test]),
      lasso=metrics(y[test],fit_lasso$fitted[test]),
      enet=metrics(y[test],fit_enet$fitted[test])
)
#selection TPR FPR
rbind(sgvc.sparse=sel_rates(fit_sgvc.sparse$sel_ind,truth,p),
      lasso=sel_rates(fit_lasso$sel_ind,truth,p),
      enet=sel_rates(fit_enet$sel_ind,truth,p)
)
