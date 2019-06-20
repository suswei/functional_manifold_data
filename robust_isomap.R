

robust_iso<- function(X){
  n <- dim(X)[1]
  D <- as.matrix(dist(X))
  A <- manif(X)
  DA <- D*A
  G <- graph.adjacency(DA, "undirected", weighted = TRUE)
  GD <- shortest.paths(G)
}


mst <- function(X) {
  D <- as.matrix(dist(X))
  G <- graph.adjacency(D, "undirected", weighted = TRUE)
  T <- minimum.spanning.tree(G)
  A <- get.adjacency(T)
  return(A)
}


# Function to obtain an extended graph
manif <- function(X) {
  D <- as.matrix(dist(X))
  A <- mst(X)
  DA <- D*A
  R <- apply(DA, 1, max)
  n <- dim(X)[1]
  d <- dim(X)[2]
  A <- matrix(0, nrow=n, ncol=n)
  L <- 100
  f <- function(l,k) {
    (1-l)*X[i,k]+l*X[j,k]
  }
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      P <- outer(1:L/(L+1), 1:d, f)
      E <- as.matrix(dist(rbind(P, X)))
      F <- E[(L+1):(L+n), 1:L]
      if (all(apply(F<R, 2, any))) A[i,j]=1
    }
  }
  A <- A + t(A)
  return(A)
}


