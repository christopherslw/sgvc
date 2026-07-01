# sparse implementation of spectral graph varying coefficient model
# group lasso chooses between a zero, a constant, and a varying effect per feature
# requires grpreg

library(grpreg)

sgvc.sparse <- function(X, y, L, train, nb) {
  p = ncol(X)
  feats = colnames(X)
  Phi = graph_basis(L, nb)
  D = Phi
  cn = paste0("phi", 1:nb)
  group = rep(1, nb)
  g = 1
  for (j in 1:p) {
    D = cbind(D, X[,j], X[,j]*Phi)
    cn = c(cn, paste0(feats[j],".c"),paste0(feats[j],".d",1:nb))
    group = c(group, g+1, rep(g+2, nb))
    g = g+2
  }
  colnames(D) = cn
  
  f = grpreg(D[train,], y[train], group=group, penalty="grLasso")
  rss = colSums((y[train]-predict(f, D[train,]))**2)
  bic = length(train)*log(rss/length(train))+f$df*log(length(train))
  jbest = which.min(bic)
  b = coef(f)[, jbest]
  names(b) = c("(Intercept)", cn)
  list(fitted=as.numeric(predict(f, D)[, jbest]), coef=b, nb=nb)
}
