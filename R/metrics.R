# evaluate different methods
metrics <- function(y, yhat){
  e = y - yhat
  R2=1-sum(e**2)/sum((y-mean(y))**2)
  RMSE=sqrt(mean(e**2))
  MAE=mean(abs(e))
  c(R2=R2,RMSE=RMSE,MAE=MAE)
}