rm(list = ls())

source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('pairwiseDistances.R')
source('full_geo_from_adj_geo.R')
source('sim_study.R')

library(fields)
library(reticulate)
library(RandPro)

rp = FALSE

# use right version of python
use_python('/anaconda3/bin/python',required=TRUE)
pyIso = import_from_path("getIsomapGdist",path='.')
scms = import_from_path("scms",path='.')

################### SIMULATION STUDY

# set up parameters

name = c("archimedean-spiral",'circle',"right-angle","sin-cos-curve")
samplesize = c(200,500)
SNR = c(20,25,30) # noise sd
reg_sampling=c(TRUE,FALSE)

# Generate the data
##### YOU NEED TO CREATE A FOLDER IN THE CURRECNT DIRRECTORY CALLED simulated_data BEFORE USING THE FUNCTION
for(mani_type in 1:length(name)){
  for(n in 1:length(samplesize)){
    for(level_noise in 1:length(SNR)){
      for(sampling in 1:length(reg_sampling)){
        generate_data(name[mani_type],samplesize[n],SNR[level_noise],reg_sampling[sampling])
      }
    }
  }
}

#Analyse the data and obtain the desired plot with do_plot and do_pdf options

for(mani_type in 1:length(name)){
  for(n in 1:length(samplesize)){
    for(level_noise in 1:length(SNR)){
      for(sampling in 1:length(reg_sampling)){
        sim_study(name[mani_type],samplesize[n],SNR[level_noise],reg_sampling[sampling],do_plot=FALSE,do_pdf=TRUE)
      }
    }
  }
}


################### END SIMULATION STUDY


# set up parameters
name = "sin-cos-curve"
samplesize = 200
noi_SNR = 24 # noise sd
scms_h = .7 # bandwidth parameter
num_neigh = 5
reg_sampling=FALSE

# graphics
par(mfrow = c(2,3))

# load true and noisy manifold data from EuclideanExamples
obj = EuclideanExamples(name, samplesize, noi_SNR, plotTrue=TRUE,reg_sampling)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)
true_geo = data.matrix(obj$true_geo)

# manifold estimation via subspace constrained mean shift (SCMS)
# TODO Susan: how to choose bandwidth?

denoised = scms$scms(data, scms_h)  #has same shape as data
if(rp==TRUE){
  # hit denoised with random projection matrix, project to same dimension
  rpmat = form_matrix(rows=2, cols=2, JLT=FALSE, eps = 0.1, projection = "gaussian")
  denoised = denoised %*% rpmat
}
plot(denoised, pch=19, xlab='', ylab='', main=paste("SCMS", 'h=', scms_h, sep=' '))



# use SKLEARN isomap
IsomapGdist_true = pyIso$getIsomapGdist(true_mani, num_neigh)
IsomapGdist_obs = pyIso$getIsomapGdist(data, num_neigh)
IsomapGdist_denoised = pyIso$getIsomapGdist(denoised, num_neigh)

# extract just upper triangular, don't include diagonal zeros
IsomapGdist_true = IsomapGdist_true[lower.tri(IsomapGdist_true, diag = FALSE)]
IsomapGdist_obs = IsomapGdist_obs[lower.tri(IsomapGdist_obs, diag = FALSE)]
IsomapGdist_denoised = IsomapGdist_denoised[lower.tri(IsomapGdist_denoised, diag = FALSE)]

# standardize geodesic distance matrices, okay to overwrite...?
# IsomapGdist_true = cov2cor(IsomapGdist_true)
# IsomapGdist_obs = cov2cor(IsomapGdist_obs)
# IsomapGdist_denoised = cov2cor(IsomapGdist_denoised)

# average Frobenius error of geodesic distance matrix
gdist_froberr_obs = mean((IsomapGdist_true - IsomapGdist_obs)^2)
gdist_froberr_denoised = mean((IsomapGdist_true - IsomapGdist_denoised)^2)
# TODO Marie: add measure assessing near-isometry

plot(IsomapGdist_true, IsomapGdist_obs, main=paste('err=', round(gdist_froberr_obs,digits=2),sep=''),
     xlab = "isomap geodesic noiseless", ylab = "isomap geodesic noisy")
abline(0,1,col="red")

plot(IsomapGdist_true,IsomapGdist_denoised, main=paste('err=', round(gdist_froberr_denoised,digits=2),sep=''),
       xlab = "isomap geodesic noiseless", ylab = "isomap geodesic scms")
abline(0,1,col="red")


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

