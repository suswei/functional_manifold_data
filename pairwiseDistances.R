## Computes the pairwise Euclidean distances between all data samples
##
## input:
##       data: NxD matrix (N samples, D features)
## output:
##	NxN matrix with the pairwise Euclidean distances


pairwiseDistances = function(data){

        num_samples = nrow(data)

        d = t(rowSums(data^2))
        d = d[rep(1, each = num_samples),]
        d = d + t(d) - 2*data%*%t(data)
        diag(d) = 0
        d = sqrt(d)

	return(d)

}
