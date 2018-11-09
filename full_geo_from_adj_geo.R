# calculate pairwise geodesic from adjacent geodesic geo(X_1,X_2),geo(X_2,X_3),etc.

full_geo <- function(adj_vec,n){
	
	
	geo_mani <- matrix(0,ncol=n,nrow=n)
	diag(geo_mani[-n,-1]) <- adj_vec
	temp <- adj_vec[-(n-1)]
	for(i in 2:(n-2)){
		
		dia_temp <- adj_vec[-(1:(i-1))]+temp
		diag(geo_mani[1:(n-i),-(1:i)]) <- dia_temp
		temp <- dia_temp[-(n-i)]
	}
	geo_mani[1,n] <- temp+adj_vec[n-1]
	return(geo_mani+t(geo_mani))
}


arc_length_spiral <- function(theta){
	
	0.5*(theta*sqrt(1+theta^2) + log(theta+sqrt(1+theta^2)))
	
}