# smooth graph basis corresponding to the nb lowest-frequency eigenvectors of the normalized Laplacian L
# eigendecomposition depends only on the graph, so it is computed once and cached
# returns an n x nb matrix
# improved for large graphs

library(RSpectra)

graph_basis <- local({
  L0 = NULL
  V = NULL
  function(L, nb) {
    n = nrow(L)
    if (is.null(L0) || !identical(L, L0) || ncol(V) < nb+1) {
      L0 <<- L
      if (inherits(L, "sparseMatrix")) {
        k = min(nb+6, n-1)
        e = eigs_sym(L, k, sigma=-1e-6) # k smallest eigenpairs
        V <<- e$vectors[, order(e$values), drop=FALSE]
      } else {
        V <<- eigen(as.matrix(L),symmetric=TRUE)$vectors[,n:1,drop=FALSE]
      }
    }
    Phi = V[, 2:(nb+1), drop=FALSE]*sqrt(n)
    colnames(Phi) = paste0("phi",1:nb)
    Phi
  }
})

