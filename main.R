source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('pairwiseDistances.R')

library(fields)
library(reticulate)

# use right version of python
use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)

# set up parameters
name = "manifold"
samplesize = 500
noi_sd <- 1.5 # noise sd
scms_h = 0.7 # bandwidth para.
# TODO Susan: why is isomap geodesic heat map so sensitive to num_neigh???
num_neigh = 11

# graphics
par(mfrow = c(2,3))

set.seed(1234)
# load true and noisy manifold data from EuclideanExamples
obj = EuclideanExamples(name, samplesize,noi_sd)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)

# manifold estimation via subspace constrained mean shift (SCMS)
# TODO Susan: how to choose bandwidth?
scms = import_from_path("scms",path='.')
denoised = scms$scms(data, scms_h)  #has same shape as data
plot(denoised, pch=19, xlab='', ylab='', main=paste("SCMS", 'h=', scms_h, sep=' '))

# use SKLEARN isomap
pyIso = import_from_path("getIsomapGdist",path='.')
IsomapGdist_true <- pyIso$getIsomapGdist(true_mani, num_neigh)
IsomapGdist_obs <- pyIso$getIsomapGdist(data, num_neigh)
IsomapGdist_denoised <- pyIso$getIsomapGdist(denoised, num_neigh)

# standardize geodesic distance matrices, okay to overwrite...?
# IsomapGdist_true = cov2cor(IsomapGdist_true)
# IsomapGdist_obs = cov2cor(IsomapGdist_obs)
# IsomapGdist_denoised = cov2cor(IsomapGdist_denoised)

# average Frobenius error of geodesic distance matrix
gdist_froberr_obs <- mean((IsomapGdist_true - IsomapGdist_obs)^2)
gdist_froberr_denoised <- mean((IsomapGdist_true - IsomapGdist_denoised)^2)
# TODO Marie: add measure assessing near-isometry

image.plot(IsomapGdist_true,main='isomap geodesic noiseless')
image.plot(IsomapGdist_obs,main=paste('isomap geodesic obs, err=',round(gdist_froberr_obs,digits=6),sep=''))
image.plot(IsomapGdist_denoised,main=paste('isomap geodesic scms, err=',round(gdist_froberr_denoised,digits=6),sep=''))


# TODO future: how to project onto
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

