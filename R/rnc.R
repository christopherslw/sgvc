# regression with network cohesion (Li et al., 2019)
# a smooth varying intercept on the graph + constant slopes per feature
#  intercept is expanded in the smooth graph basis Phi
# rewritten using logic from (Li et al., 2019) paper and netcoh R package

rnc <- function(X, y, L, train=NULL,lambdas=NULL,nfolds=5, gamma=0.01) {
  if (is.null(train)){
    train = 1:length(y)
  }
  if (is.null(lambdas)){
    lambdas=10^seq(-4,4,length.out=20)
  }
  X = as.matrix(X)
  L = as.matrix(L)
  n = nrow(X)
  p = ncol(X)
  Lg = L + gamma*diag(n)
  
  fit_one = function(obs, lambda) {
    Pv = numeric(n); Pv[obs] = 1
    PX = X*Pv
    Py = y*Pv
    C = lambda*Lg
    diag(C) = diag(C)+Pv 
    Z = solve(C, cbind(PX, Py))
    Zx = Z[, 1:p, drop=FALSE]; zy = Z[, p+1]
    S = crossprod(X[obs,,drop=FALSE]) - crossprod(PX, Zx)   # Schur complement
    rhs = crossprod(X[obs,,drop=FALSE], y[obs]) - crossprod(PX, zy)
    beta = solve(S, rhs)
    alpha = zy-Zx%*%beta
    list(alpha=as.numeric(alpha), beta=as.numeric(beta))
  }
  
  # choose lambda by k fold CV 
  fold = sample(rep(1:nfolds, length.out=length(train)))
  cve = matrix(NA, nfolds, length(lambdas))
  for (k in 1:nfolds) {
    obs = train[fold != k]
    hold = train[fold == k]
    for (l in seq_along(lambdas)) {
      ft = fit_one(obs, lambdas[l])
      yh = ft$alpha+X%*%ft$beta
      cve[k, l] = mean((y[hold]-yh[hold])^2)
    }
  }
  lambda = lambdas[which.min(colMeans(cve))]
  
  ft = fit_one(train, lambda)
  beta_hat = ft$beta
  names(beta_hat) = colnames(X)
  fitted = as.numeric(ft$alpha+X%*%ft$beta)
  return(list(fitted=fitted,alpha_hat=ft$alpha, beta_hat=beta_hat,
              lambda=lambda, lambdas=lambdas, cve=colMeans(cve)))
}