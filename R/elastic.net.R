# requires glmnet
library(glmnet)

elastic.net <- function(X, y, train=NULL){
  if(is.null(train)){
    train = 1:length(y)
  }
  fit = cv.glmnet(X[train,], y[train], alpha=0.5)
  fitted= as.numeric(predict(fit, X, s="lambda.min"))
  beta_hat = as.numeric(coef(fit, s="lambda.min"))[-1]
  names(beta_hat) = colnames(X)
  sel_ind = which(beta_hat != 0)
  nsel = length(sel_ind)
  
  return(list(fitted=fitted, beta_hat=beta_hat, nsel=nsel, sel_ind=sel_ind))
}
