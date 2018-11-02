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
name = "spiral"
samplesize = 2000
num_neigh = 5


par(mfrow = c(1,2))

# extract irregularly sampled true manifold
#set.seed(1235)
noi_sd = 0.1
obj = EuclideanExamples(name, samplesize, noi_sd, plotTrue = FALSE)
true_mani  = data.matrix(obj$true_mani)
plot(true_mani, pch=1, xlab='', ylab='')

# calculate geodesic distances of true manifold, check for collapsing
IsomapGdist_truemani = pyIso$getIsomapGdist(true_mani, num_neigh)
collapseIndex = which(IsomapGdist_truemani == 0, arr.ind = TRUE)
collapseIndex= collapseIndex[ collapseIndex[,1]!=collapseIndex[,2], ]
print(length(collapseIndex))
if(length(collapseIndex)>0)
{
  k=2
  ind1 = collapseIndex[k,1]
  print(ind1)
  ind2 = collapseIndex[k,2]
  print(ind2)
  points(true_mani[ind1,1],true_mani[ind1,2],col="blue")
  points(true_mani[ind2,1],true_mani[ind2,2],col="blue")
}

# histogram of nondiagonal pairwise geodesic differences
IsomapGdist_truemani = IsomapGdist_truemani[lower.tri(IsomapGdist_truemani, diag = FALSE)]
hist(IsomapGdist_truemani)




