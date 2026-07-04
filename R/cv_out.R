# k-fold cross validation  
# method: any function of the training indices that returns a list with a fitted vector of length n
# returns the out of sample predictions, i.e. each node's value is taken from the fit of the fold in which that node was held out.

cv_out <- function(method, fold){
  n = length(fold)
  yhat = numeric(n)
  for (f in sort(unique(fold))) {
    tr = which(fold != f)
    yhat[fold == f] = method(tr)$fitted[fold == f]
  }
  yhat
}
