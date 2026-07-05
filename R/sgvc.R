# sparse implementation of spectral graph varying coefficient model
# group lasso chooses between a zero, a constant, and a varying effect per feature
# requires grpreg
library(grpreg)

sgvc <- function(X, y, L, train=NULL, nb=10, screen=500, nfolds=5) {
  if (is.null(train)){
    train = 1:length(y)
  }
  n = nrow(X)
  p = ncol(X)
  Phi = graph_basis(L, nb)

  # screen group score of each feature against centered y
  r = numeric(n)
  r[train] = y[train]-mean(y[train])
  score = rowSums(crossprod(X, cbind(r, Phi * r))^2)
  keep = sort(order(score, decreasing=TRUE)[1:min(screen, p)])
  k = length(keep)

  D = matrix(0, n, nb+k*(nb+1))
  group = integer(ncol(D))
  D[,1:nb] = Phi
  group[1:nb] = 1
  col = nb
  for (i in 1:k) {
    D[, col+1] = X[,keep[i]]
    D[, col+1 + (1:nb)] = X[,keep[i]]*Phi
    group[col+1] = 2*i
    group[col+1+(1:nb)] = 2*i+1
    col = col+nb+1
  }

  # fit path on training nodes & choose lambda
  f = grpreg(D[train,], y[train], group=group, penalty="grLasso")
  lam = f$lambda
  fold = sample(rep(1:nfolds, length.out=length(train)))
  mse = matrix(NA, nfolds, length(lam))
  for (v in 1:nfolds) {
    tr = train[fold!=v]
    te = train[fold==v]
    g = grpreg(D[tr,], y[tr], group=group, penalty="grLasso", lambda=lam)
    pred = predict(g, D[te,,drop=FALSE])
    mse[v, 1:ncol(pred)]=colMeans((y[te]-pred)^2)
  }
  best = which.min(colMeans(mse, na.rm=TRUE))

  b = coef(f)[, best]
  fitted = as.numeric(predict(f, D, which=best))
  B = matrix(b[-(1:(nb + 1))], nb + 1, k)
  sel_ind = keep[colSums(B != 0) > 0]
  names(sel_ind) = colnames(X)[sel_ind]

  return(list(fitted=fitted, coef=b, nb=nb, lambda=lam[best], nsel=length(sel_ind),
              sel_ind=sel_ind, screened=keep, f=f))
}
