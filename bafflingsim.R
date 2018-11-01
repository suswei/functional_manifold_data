# TODO: this simulation exhibits collapsing ISOMAP geodesic distance, what is this due to?
# Seems like it's pretty easy to find other cases of collapse in the Euclidean examples for certain sd_noise and scms_h
# related to short circuit?

rm(list = ls())

source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('pairwiseDistances.R')

library(fields)
library(reticulate)
library(RandPro)

### blaha blah blah

rp = TRUE
# use right version of python
use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)
pyIso = import_from_path("getIsomapGdist",path='.')

# set up parameters
name = "archimedean-spiral"
samplesize = 200
scms_h = .5 # bandwidth parameter
num_neigh = 5

par(mfrow = c(3,2))

# noise at 0.2
noi_sd = 0.2 # noise sd
obj = EuclideanExamples(name, samplesize, noi_sd, plotTrue = FALSE)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)
scms = import_from_path("scms",path='.')
denoised02 = scms$scms(data, scms_h)

if(rp==TRUE){
  # hit denoised with random projection matrix, project to same dimension
  rpmat = form_matrix(rows=2, cols=2, JLT=FALSE, eps = 0.1, projection = "gaussian")
  denoised02 = denoised02 %*% rpmat
}

# check if scms collapses many points, vector of unique pairwise distances in denoised02
dists02 = dist(denoised02, upper=TRUE)
hist(dists02, main="pariwise distances in denoised02")
# apply isomap geodesic distance subroutine
IsomapGdist_denoised02 = pyIso$getIsomapGdist(denoised02, num_neigh)
IsomapGdist_denoised02 = IsomapGdist_denoised02[lower.tri(IsomapGdist_denoised02, diag = FALSE)]

######################
set.seed(1011) # have seen collapsing with seed 1011 and 1234
# noise at 0.3
noi_sd = 0.3 # noise sd
obj = EuclideanExamples(name, samplesize,noi_sd, plotTrue = FALSE)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)
scms = import_from_path("scms",path='.')
denoised03 = scms$scms(data, scms_h)

if(rp==TRUE){
  # hit denoised with random projection matrix, project to same dimension
  rpmat = form_matrix(rows=2, cols=2, JLT=FALSE, eps = 0.1, projection = "gaussian")
  denoised03 = denoised03 %*% rpmat
}

# check if scms collapses many points, vector of unique pairwise distances in denoised02
dists03 = dist(denoised03, upper=TRUE)
hist(dists03, main="pariwise distances in denoised03")
# apply isomap geodesic distance subroutine
IsomapGdist_denoised03 = pyIso$getIsomapGdist(denoised03, num_neigh)
IsomapGdist_denoised03 = IsomapGdist_denoised03[lower.tri(IsomapGdist_denoised03, diag = FALSE)]

######################
plot(denoised02, pch=19, xlab='', ylab='', main=paste("SCMS", 'h=', scms_h, sep=' '))
points(denoised03, pch=19, col = "red")

plot(IsomapGdist_denoised02,IsomapGdist_denoised03,xlab='IsomapGdist_denoised02', ylab = 'IsomapGdist_denoised03')
abline(0,1,col="red")

hist(IsomapGdist_denoised02)
hist(IsomapGdist_denoised03)
