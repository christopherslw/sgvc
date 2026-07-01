# requires glmnet
library(glmnet)
lasso <- function(X, y, train)
{
  fit = cv.glmnet(X[train,], y[train], alpha=1)
  list(fitted=as.numeric(predict(fit, X, s="lambda.min")))
}
