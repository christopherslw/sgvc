# stochastic block model graph with K equal communities
sbm_graph <- function(n, K=4, pin=0.12, pout=NULL) {
  if (is.null(pout)) pout <- 4 / n
  sizes <- rep(floor(n / K), K)
  sizes[K] <- n - sum(sizes[-K])
  block <- rep(seq_len(K), times=sizes)
  ii <- integer(0)
  jj <- integer(0)
  for (a in seq_len(n - 1)) {
    nb_block <- block[(a + 1):n]
    pr  <- ifelse(nb_block == block[a], pin, pout)
    hit <- which(runif(length(pr)) < pr)
    if (length(hit)) {
      ii <- c(ii, rep(a, length(hit)))
      jj <- c(jj, a + hit)
    }
  }
  A <- sparseMatrix(i=ii, j=jj, x=1, dims=c(n, n))
  list(W=drop0(1 * ((A + t(A)) > 0)), block=block)
}

# community coefficients under the SBM
# constant within each block
comm_coef <- function(block, s, level_sd=1) {
  K <- max(block)
  B <- matrix(0, length(block), s)
  for (j in seq_len(s)) B[, j] <- rnorm(K, 0, level_sd)[block]
  B
}


# 2-D lattice graph
lattice_graph <- function(side) {
  n <- side * side
  idx <- function(r, c) (r - 1) * side + c
  ii <- integer(0)
  jj <- integer(0)
  for (r in seq_len(side)) for (c in seq_len(side)) {
    v <- idx(r, c)
    if (c < side) {
      ii <- c(ii, v)
      jj <- c(jj, idx(r, c + 1))
    }
    if (r < side) {
      ii <- c(ii, v)
      jj <- c(jj, idx(r + 1, c))
    }
  }
  A <- sparseMatrix(i=ii, j=jj, x=1, dims=c(n, n))
  coords <- cbind(rep(seq_len(side),each=side),rep(seq_len(side),times=side))/side
  list(W=drop0(1 * ((A + t(A)) > 0)), coords=coords, n=n)
}

#smooth functions of 2-D coordinates
smooth_lib <- list(
  function(c)  0.8 + 1.2 * sin(pi * c[, 1]) * sin(pi * c[, 2]),
  function(c) -0.5 + 1.5 * cos(pi * c[, 1]),
  function(c)  1.0 * (c[, 1] - 0.5) + 1.0 * (c[, 2] - 0.5),
  function(c)  1.5 * exp(-((c[, 1] - .5)^2 + (c[, 2] - .5)^2) / 0.15) - 0.5,
  function(c)  0.6 + 1.0 * cos(2 * pi * c[, 2]),
  function(c)  1.2 * (c[, 1]^2 - c[, 2]^2)
)


# generate simulation data
gen_data <- function(n, p, s, snr=5, graph=c("sbm", "lattice"),
                     rho_x=0, tau=2.5, vary=0) {
  graph <- match.arg(graph)
  coords <- NULL
  block <- NULL
  if (graph == "lattice") {
    lat <- lattice_graph(round(sqrt(n)))
    W <- lat$W
    coords <- lat$coords
    n <- lat$n
  } else {
    sg <- sbm_graph(n)
    W <- sg$W
    block <- sg$block
  }
  Z <- matrix(rnorm(n * p), n, p) # covariates
  if (rho_x > 0) {
    X <- Z
    for (j in 2:p){
      X[, j] <- rho_x*X[, j-1] + sqrt(1-rho_x^2)*Z[, j]
    }
  } else X <- Z
  X <- scale(X)
  active <- sort(sample(p, s))
  Bstar <- matrix(0, n, p)
  alpha_star <- numeric(n)
  if (!is.null(coords)) {
    for (ti in seq_along(active)) {
      col <- smooth_lib[[((ti - 1) %% length(smooth_lib)) + 1]](coords)
      Bstar[, active[ti]] <- col / sd(col)
    }
  }
  else {
    b <- rnorm(s)
    b <- b * sqrt(s / sum(b^2)) 
    Bstar[, active] <- matrix(b, n, s, byrow=TRUE)
    if (vary > 0)
      Bstar[, active] <- Bstar[, active] + vary * comm_coef(block, s)
    a <- as.numeric(comm_coef(block, 1))
    alpha_star <- tau * a / sd(a)
  }
  fstar <- alpha_star + rowSums(X * Bstar)
  sigma <- sd(fstar)/sqrt(snr)
  y <- fstar + rnorm(n, 0, sigma)
  list(X=X, y=y, W=W, coords=coords, block=block, Bstar=Bstar, alpha_star=alpha_star,
       fstar=fstar, active=active, sigma=sigma, n=n, p=p, s=s)
}