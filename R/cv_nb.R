# choose the spectral basis length nb by k-fold cross validation
# candidate lengths are tried in increasing order and out of sample error is formed with cv_out on the shared folds.
# since the error is smooth in nb, search stops once it has failed to improve for 'patience'

cv_nb <- function(fit, grid, fold, y, patience=3)
{
  best = Inf
  best_nb = grid[1]
  worse = 0
  for (nb in grid) {
    err = mean((y - cv_out(function(tr) fit(tr, nb), fold)) ** 2)
    if (err < best) {
      best = err
      best_nb = nb
      worse = 0
    } else {
      worse = worse + 1
      if (worse >= patience) break
    }
  }
  best_nb
}
