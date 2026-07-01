# regression with network cohesion (Li, Levina & Zhu, 2019)
# a smooth varying intercept on the graph + constant slopes per feature
#  intercept is expanded in the smooth graph basis Phi
# incomplete
rnc <- function(X, y, Phi, train)
{
  M = cbind(1, Phi, X)
  b = solve(crossprod(M[train,]) + 1e-4*diag(ncol(M)), crossprod(M[train,], y[train]))
  list(fitted=as.numeric(M %*% b))
}
