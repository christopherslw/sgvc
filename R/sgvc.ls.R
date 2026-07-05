# fit the spectral graph varying coefficient model
# returns the fitted values and the coefficient functions
# fit obtained by qr rather than the normal equations, which is more robust when an eigenvector is nearly colinear with a covariate

sgvc.ls <- function(X, y, L, train=NULL, nb=NULL, nb_max=20, eps=1e-4) {
  if (is.null(train)) train = 1:length(y)
  p = ncol(X)
  press = NULL
  
  if (is.null(nb)) {
    K = min(nb_max, floor((length(train)-2-p)/(p+1)))
    if (K < 1){
      stop("training set too small")
    }
    Phi = graph_basis(L, K)
    M = cbind(1, X)
    for (l in 1:K){
      M = cbind(M,Phi[,l],X*Phi[,l])
    }
    d = ncol(M)
    qa = qr(rbind(M[train,], sqrt(eps)*diag(d)))
    z  = qr.qty(qa, c(y[train], numeric(d)))[1:d]
    Q  = qr.Q(qa)[1:length(train),,drop=FALSE]
    Qz = Q*rep(z, each=length(train))
    Q2 = Q^2
    
    press = sapply(1:K, function(m) {
      k = 1 + p + m*(1+p)
      fitted = rowSums(Qz[, 1:k, drop=FALSE])
      leverage = rowSums(Q2[, 1:k, drop=FALSE])
      mean(((y[train]-fitted)/(1-leverage))^2)
    })
    nb = which.min(press)
  }
  
  Phi = graph_basis(L, nb)
  D = Phi
  for (j in 1:p) D = cbind(D, X[,j], X[,j]*Phi)
  M = cbind(1, D)
  d = ncol(M)
  Ma = rbind(M[train,],sqrt(eps)*diag(d))
  b = qr.coef(qr(Ma), c(y[train], numeric(d)))
  list(fitted=as.numeric(M%*%b), coef=b, nb=nb, press=press)
}