# choose the spectral basis length nb for the least squares fit by k-fold cross validation
# leverages existing basis and only compute one qr factorization per fold 

sgvc.ls.cv <- function(X, y, L, fold, grid, patience=3)
{
  n = nrow(X)
  p = ncol(X)
  Phi = graph_basis(L, max(grid))
  M = cbind(1, X)
  for (l in 1:max(grid))
    M = cbind(M, Phi[,l], X*Phi[,l])
  bw = p+1
  d = ncol(M)
  pred = matrix(0, n, max(grid))
  for (f in sort(unique(fold))) {
    tr = which(fold != f)
    ho = which(fold==f)
    Ma = rbind(M[tr,], sqrt(1e-4)*diag(d))
    qro = qr(Ma)   # one factorization per fold reused across nb
    R = qr.R(qro)
    z = qr.qty(qro, c(y[tr], numeric(d)))[1:d]
    for (nb in grid) {
      nc = bw * (nb + 1)
      b = backsolve(R[1:nc, 1:nc], z[1:nc])
      pred[ho, nb] = M[ho, 1:nc]%*%b
    }
  }
  err = sapply(grid, function(nb) mean((y-pred[, nb])**2))
  best = Inf
  best_nb = grid[1]
  worse = 0
  for (i in seq_along(grid)) {
    if (err[i] < best) {
      best = err[i]
      best_nb = grid[i]
      worse = 0
    } else {
      worse = worse+1
      if (worse >= patience) break
    }
  }
  best_nb
}
