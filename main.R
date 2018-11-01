source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('pairwiseDistances.R')

library(fields)
library(reticulate)
use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)
scms = import_from_path("scms",path='.')
par(mfrow=c(1,3))

# load noisy manifold data from EuclideanExamples
name = "manifold"
samplesize = 1000
obj <- EuclideanExamples(name, samplesize)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)

#set up parameters
name = "manifold"
samplesize = 500
noi_sd <- 1.5 #noise sd
scms_h = 1 #bandwidth para.
dim_mani <- 1
num_neigh <- 10

# manifold estimation via local covariance flow
# x0 = as.numeric(obj$data[500,])
# sigma = 1
# kernel = scms$make_isotropic_gaussian_kernel(sigma)
# ms_iterations = 1000
# for (j in 1:ms_iterations){
#   x0shifted = x0 + scms$mean_shift_update(x0, data, kernel)
#   x0 = x0shifted
# }
# shift_type = "ms"
# kernel = "double-peak"
# eig_h = 0.9
# path = eigenvectorFlow(obj, x0shifted, shift_type, kernel, eig_h, FALSE)




par(mfrow=c(3,3))

# load noisy manifold data from EuclideanExamples
obj <- EuclideanExamples(name, samplesize,noi_sd)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)


# manifold estimation via subspace constrained mean shift (SCMS)

denoised = scms$scms(data, scms_h)  #has same shape as data
plot(denoised, pch=19, xlab='', ylab='', main=paste("SCMS", 'h=', scms_h, sep=' '))


#Use Isomap to obtain manifold representation [[1]] and geodesic distance [[2]]
Iso_true <- Isomap(true_mani, dim_mani, num_neigh)
Iso_obs_data <- Isomap(data, dim_mani, num_neigh)
Iso_denoised_data <- Isomap(denoised, dim_mani, num_neigh)

# standardize geodesic distance matrices


rel_err_obs_data <- sum((Iso_true[[2]]-Iso_obs_data[[2]])^2)/sum(Iso_true[[2]]^2)
rel_err_denoised_data <- sum((Iso_true[[2]]-Iso_denoised_data[[2]])^2)/sum(Iso_true[[2]]^2)

image.plot(Iso_true[[2]],main='true geodesic')
image.plot(Iso_obs_data[[2]],main=paste('obs. geo., rel. err=',round(rel_err_obs_data,digits=6),sep=''))
image.plot(Iso_denoised_data[[2]],main=paste('denoi. geo, rel. err=',round(rel_err_denoised_data,digits=6),sep=''))

if(dim_mani==2){
	Iso_comp_true <- Iso_true[[1]][[1]]
	plot(Iso_comp_true[,1],Iso_comp_true[,2],main='mani comp true')
	Iso_comp_obs <- Iso_obs_data[[1]][[1]]
	plot(Iso_comp_obs[,1],Iso_comp_obs[,2],main='mani comp obs.')
	Iso_comp_den <- Iso_denoised_data[[1]][[1]]
	plot(Iso_comp_den[,1],Iso_comp_den[,2],main='mani comp denoi.')
} else if (dim_mani==1){
	Iso_comp_true <- Iso_true[[1]][[1]]
	plot(1:samplesize,Iso_comp_true[,1],main='mani comp true')
	Iso_comp_obs <- Iso_obs_data[[1]][[1]]
	plot(1:samplesize,Iso_comp_obs[,1],main='mani comp obs.')
	Iso_comp_den <- Iso_denoised_data[[1]][[1]]
	plot(1:samplesize,Iso_comp_den[,1],main='mani comp denoi.')
}

