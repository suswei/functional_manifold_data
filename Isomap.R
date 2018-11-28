## Computes the Isomap embedding by Tenenbaum, de Silva and Langford (2000).
## It is based on the Matlab implementation by Tenenbaum.
## This version uses Floyd's Algorithm to compute the shortest path tree.
##
## A modified version of the original Isomap algorithm is included. It respects nearest and farthest neighbours.
## To estimate the intrinsic dimension of the data, the function can plot the residuals between the high and the 
## low dimensional data for a given range of dimensions.
##
## input: 
##	data: 		NxD matrix (N samples, D features)
##	dims:		dimension of the target space (default 2)
##	k:		number of neighbours (default 5 or number of samples)
##	mod:		use modified Isomap algorithm (default FALSE)
##	plotResiduals:	show a plot with the residuals between the high and the low dimensional data (default FALSE)
##	verbose:	show a summary of the embedding procedure at the end (default TRUE)	 
##
## output: 
##	Nxdim matrix (N samples, dim features) with the reduced input data (list of several matrices if more than one dimension specified)


Isomap = function(data, dims=2, k, mod=FALSE){
#Isomap = function(data,  k, mod=FALSE, plotResiduals=FALSE, verbose=FALSE){
# catch missing or bad input
	if(missing(data))
      stop("data argument missing")
  else
    if(!is.matrix(data))
        stop("invalid argument: data argument is required to be a N x D matrix (N samples, D features)")
    if(!all(is.numeric(dims)) | any(dims < 1))
        stop("invalid argument: target dimension is required to be a (vector of) positive integer value(s)")
	if(missing(k))
		k = min(nrow(data),5)
	else{
    if(k >= nrow(data))
      stop("invalid argument: more neighbours than samples")
    if(!is.numeric(k) | k <= 1)
      stop("invalid argument: neighbour parameter is required to be an integer value >= 2")
  }

  k = round(k)
  dims = round(dims)
        
  num_samples = nrow(data)
  num_features = ncol(data)

  ## compute pairwise distances
	d = pairwiseDistances(data)
	
	## build graph of shortest paths between two neighbours (using Floyd's Algorithm)
	## modified Isomap: Use k/2 nearest neighbours and k/2 farthest neighbours
	if(mod == TRUE){
		first_k = round(k/2)
		last_k = k - first_k
	}
	## original Isomap: Use k nearest neighbours only
	else{
		first_k = k
		last_k = 0
	}
	sort_idx = apply(d,2,order)

	## set weights to "Infinite" for not neighboured points
	for(i in 1:num_samples)
		d[i,sort_idx[(first_k+2):(num_samples-last_k),i]] = Inf

	## ensure that graph matrix is symmetric
	d = pmin(d,t(d))

	## Floyd's Algorithm
	for(i in 1:num_samples){
		d_col_i = t(d[,i])
		d_col_i = t(d_col_i[rep(1, each = num_samples),])
		d_row_i = t(d[i,])
		d_row_i = d_row_i[rep(1, each = num_samples),]
		d = pmin(d, d_col_i + d_row_i)
	}
        
	## determine all connected components of the graph
	num_connections = rowSums(!(is.infinite(d)))
	first_connections = apply(is.infinite(d),1,which.min)
	components = unique(first_connections)
	num_components = length(components)

	## reduce dimension of all components seperately and merge them to a single dataset
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
	
	rep <- list(all_Y,d)
	names(rep) <- c("bla","geod")

	return(rep)

}
