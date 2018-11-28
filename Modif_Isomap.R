## Computes the Isomap embedding by Tenenbaum, de Silva and Langford (2000).
## It is based on the Matlab implementation by Tenenbaum.
## This version uses Floyd's Algorithm to compute the shortest path tree.
##
## A modified version of the original Isomap algorithm is included. It respects nearest and farthest neighbours.
## To estimate the intrinsic dimension of the data, the function can plot the residuals between the high and the 
## low dimensional data for a given range of dimensions.
##
## input: 
##	d: 		NxN matrix (N samples) of geo pairwise distances
##	dims:		dimension of the target space (default 2)
##
## output: 
##	Nxdim matrix (N samples, dim features) with the reduced input data (list of several matrices if more than one dimension specified)


Modif_Isomap = function(d, dims=2){
  
  dims = round(dims)
  num_samples = nrow(d)

  ## determine all connected components of the graph
  num_connections = rowSums(!(is.infinite(d)))
  first_connections = apply(is.infinite(d),1,which.min)
  components = unique(first_connections)
  num_components = length(components)
 
  ## reduce dimension of all components seperately and merge them to a single dataset
  #message("Computing low dimensional embedding ... ", appendLF = FALSE)
  all_Y = list(NULL)
  residuals = rep(0,length(dims))
  i = 1
  for(dim in dims){
    Y = matrix(0,num_samples,dim)
    for(c in 1:num_components){
      
      ## only use the distances of points in the connected component
      comp_indices = which(first_connections == components[c])
      D = d[comp_indices,comp_indices]
      N = length(comp_indices)
      
      ## convert graph matrix to inner products
      T = -0.5 * ( D^2 - rowSums(D^2) %*% t(rep(1/N,N)) - rep(1,N) %*% t(rowSums(D^2))/N + sum(D^2)/(N^2))
      
      ## catch possible errors caused by too small connected components
      if(dim(T)[1] < dim)
        stop("Problem occured while computing embedding: Connected component ", c, " consists of too few samples (", dim(T)[1], ") and cannot be reduced to ", dim, " dimensions. Try Isomap with a neighbourhood parameter higher than ", k, " or with dimension parameters lower than ", dim, ".")
      
      ## embedding Y is given by the eigenvectors of T belonging to (some of) its biggest eigenvalues
      eig_T = eigen(T)
      sweep(eig_T$vectors, 2, sqrt(colSums(eig_T$vectors^2)), "/")	#normalize eingenvectors
      Y[comp_indices,] = Re(eig_T$vectors[,1:dim] * (matrix(1,N,1) %*%  sqrt(as.complex(eig_T$values[1:dim]))))
    }
   
    all_Y[[i]] = Y
    i = i+1
  }

  return(all_Y)
  
}
