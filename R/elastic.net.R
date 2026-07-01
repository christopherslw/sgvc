# requires glmnet
library(glmnet)

elastic.net <- function(X, y, train)
{
  fit = cv.glmnet(X[train,], y[train], alpha=0.5)
  list(fitted=as.numeric(predict(fit, X, s="lambda.min")))
}
