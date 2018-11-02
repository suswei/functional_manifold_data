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

# use right version of python
use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)
pyIso = import_from_path("getIsomapGdist",path='.')

# set up parameters
name = "sin-curve"
samplesize = 200
scms_h = 2 # bandwidth parameter
num_neigh = 5

par(mfrow = c(3,2))

set.seed(1011)
# noise at 0.2
noi_sd = 0.2 # noise sd
obj = EuclideanExamples(name, samplesize, noi_sd, plotTrue = FALSE)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)
scms = import_from_path("scms",path='.')
denoised02 = scms$scms(data, scms_h)


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
obj = EuclideanExamples(name, samplesize, noi_sd, plotTrue = FALSE)
data = data.matrix(obj$data)
true_mani  = data.matrix(obj$true_mani)
scms = import_from_path("scms", path='.')
denoised03 = scms$scms(data, scms_h)

# check if scms collapses many points, vector of unique pairwise distances in denoised02
dists03 = dist(denoised03, upper=TRUE)
hist(dists03, main="pariwise distances in denoised03")
# apply isomap geodesic distance subroutine
IsomapGdist_denoised03 = pyIso$getIsomapGdist(denoised03, num_neigh)
collapseIndex = which(IsomapGdist_denoised03 == 0, arr.ind = TRUE)

temp = lower.tri(IsomapGdist_denoised03, diag = FALSE)
IsomapGdist_denoised03 = IsomapGdist_denoised03[temp]

######################
plot(denoised02, pch=1, xlab='', ylab='', main=paste("SCMS", 'h=', scms_h, sep=' '))
points(denoised03, pch=1, col = "red")
k=2
ind1 = collapseIndex[k,1]
print(ind1)
ind2 = collapseIndex[k,2]
print(ind2)
points(denoised03[ind1,1],denoised03[ind1,2],col="blue")
points(denoised03[ind2,1],denoised03[ind2,2],col="blue")


plot(IsomapGdist_denoised02,IsomapGdist_denoised03,xlab='IsomapGdist_denoised02', ylab = 'IsomapGdist_denoised03')
abline(0,1,col="red")

hist(IsomapGdist_denoised02)
hist(IsomapGdist_denoised03)
