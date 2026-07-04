# network lasso (Hallac, Leskovec & Boyd, 2015)
# graph fused penalty  lambda * sum_{(j,k) in E} w_jk ||beta_j - beta_k||
# solved by ADMM with edge splitting
# *Assumes a connected graph*

network.lasso <- function(X, y, A, train, lambdas=2**seq(-3, 5), rho=1, iters=200){
  if(is.null(train)){
    train = 1:length(y)
  }
  X = cbind(1, X) # varying intercept plus the slopes
  n = nrow(X)
  p = ncol(X)
  edge = which(A != 0 & upper.tri(A), arr.ind=TRUE)
  j = edge[,1]
  k = edge[,2]
  w = A[edge]
  m = length(j)
  deg = tabulate(c(j, k), n)
  obs = rep(0, n)
  obs[train] = 1
  xtx = rowSums(X**2)
  xy = X*(obs*y)

  admm <- function(lambda)
  {
    beta = matrix(0, n, p)
    a = matrix(0, m, p)
    b = matrix(0, m, p)
    pj = matrix(0, m, p)
    qk = matrix(0, m, p)
    for (it in 1:iters) {
      agg = rowsum(rbind(a - pj, b - qk), c(j, k))
      S = matrix(0, n, p)
      S[as.integer(rownames(agg)),] = agg
      rhs = xy + rho*S
      c0 = rho*deg
      beta = rhs / c0
      num = rowSums(X[train,,drop=FALSE]*rhs[train,,drop=FALSE])
      beta[train,] = rhs[train,,drop=FALSE] / c0[train]-(num/(c0[train]*(c0[train]+xtx[train])))*X[train,,drop=FALSE]
      sj = beta[j,,drop=FALSE] + pj
      sk = beta[k,,drop=FALSE] + qk
      diff = sj - sk
      mid = (sj + sk) / 2
      nrm = sqrt(rowSums(diff**2))
      shrink = pmax(0, 1 - (2*lambda*w/rho)/nrm)
      delta = shrink * diff
      a = mid + delta / 2
      b = mid - delta / 2
      pj = pj + beta[j,,drop=FALSE] - a
      qk = qk + beta[k,,drop=FALSE] - b
    }
    beta
  }

  bic_min = Inf
  fitted = NULL
  for (lambda in lambdas) {
    beta = admm(lambda)
    yhat = rowSums(X * beta)
    rss = sum((y[train] - yhat[train]) ** 2)
    df = nrow(unique(round(beta, 2))) * p
    bic = length(train)*log(rss / length(train)) + df*log(length(train))
    if (bic < bic_min) {
      bic_min = bic
      fitted = yhat
    }
  }
  return(list(fitted=fitted))
}
