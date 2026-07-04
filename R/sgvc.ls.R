# fit the spectral graph varying coefficient model
# returns the fitted values and the coefficient functions
# fit obtained by qr rather than the normal equations, which is more robust when an eigenvector is nearly colinear with a covariate
# need to fix

sgvc.ls <- function(X, y, L, train, nb){
  if(is.null(train)){
    train = 1:length(y)
  }
  p = ncol(X)
  Phi = graph_basis(L, nb)
  D = Phi
  for (j in 1:p)
    D = cbind(D, X[,j], X[,j]*Phi)
  M = cbind(1, D)
  d = ncol(M)
  Ma = rbind(M[train,], sqrt(1e-4)*diag(d))
  b = qr.coef(qr(Ma), c(y[train], numeric(d)))
  list(fitted=as.numeric(M%*%b), coef=b, nb=nb)
}