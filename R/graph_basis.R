# smooth graph basis corresponding to the nb lowest-frequency eigenvectors of the normalized Laplacian L
# eigendecomposition depends only on the graph, so it is computed once and cached
# returns an n x nb matrix
# (SLOW FOR LARGE GRAPHS - may need something like a Chebyshev polynomial approximation)

graph_basis <- local({
  L0 = NULL
  eig = NULL
  function(L, nb) {
    if (is.null(L0) || !identical(L, L0)) {
      L0 <<- L
      eig <<- eigen(L, symmetric=TRUE)  #cache
    }
    n = nrow(L)
    Phi = eig$vectors[, (n - 1):(n - nb), drop=FALSE] * sqrt(n)
    colnames(Phi) = paste0("phi", 1:nb)
    Phi
  }
})
