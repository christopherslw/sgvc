# ordinary least squares
ols <- function(X, y, train=NULL) {
  if (is.null(train)) {
    train = 1:length(y)
  }
  fit = lm.fit(cbind(1, X[train,,drop=FALSE]), y[train])
  b = fit$coefficients
  b[is.na(b)] = 0
  beta_hat = b[-1]
  names(beta_hat) = colnames(X)
  fitted = as.numeric(b[1] + X %*% beta_hat)
  sel_ind = which(beta_hat != 0) 
  nsel = length(sel_ind)
  list(fitted=fitted, beta_hat=beta_hat, intercept=b[1],
       nsel=nsel, sel_ind=sel_ind)
}
