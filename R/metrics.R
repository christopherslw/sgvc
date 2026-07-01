# evaluate different methods
metrics <- function(yhat){
  e = y - yhat
  c(R2=1-sum(e**2)/sum((y-mean(y))**2), 
    RMSE=sqrt(mean(e ** 2)), 
    MAE=mean(abs(e))
    )
}