<<<<<<< HEAD
=======
setwd("/Users/marie_uqam/Dropbox/Marie-Moi/Susan_project/manifold learning/ridge_density_manifold_learning")

>>>>>>> MH_code
source('EuclideanExamples.R')
source('functions.R')

<<<<<<< HEAD

=======
>>>>>>> MH_code
library(reticulate)
use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)
scms = import_from_path("scms",path='.')
<<<<<<< HEAD
par(mfrow=c(1,3))

# load noisy manifold data from EuclideanExamples
name = "manifold"
samplesize = 1000
obj <- EuclideanExamples(name, samplesize)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)

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


# manifold estimation via subspace constrained mean shift (SCMS)
scms_h = 3

denoised = scms$scms(data, scms_h)  #has same shape as data
plot(denoised, pch=19, xlab='', ylab='', main=paste("SCMS", 'h=', scms_h, sep=' '))
