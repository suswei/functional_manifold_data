rm(list = ls())

# setwd("/Users/marie_uqam/Dropbox/Marie-Moi/Susan_project/manifold learning/Geo_dist_calculation")

source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('Modif_Isomap.R')
source('pairwiseDistances.R')
source('full_geo_from_adj_geo.R')
source('smooth_FD_bspline.R')
source('sim_functional_data.R')
source('Geo_estimation.R')
source('sim_Euclidean_datagit.R')
source('assess_goodness_estimation.R')

library(fields)
library(reticulate)
library(fda)
library(matlabr)

# use right version of python
# use_python('/anaconda3/bin/python',required=TRUE)
# use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)
pyIso = import_from_path("getIsomapGdist",path='.')
py_min_neigh = import_from_path("get_min_num_neighbors",path='.')
scms = import_from_path("scms",path='.')



#### Analysis of Euclidean data

FD_true =FALSE

# set up parameters

name = "archimedean-spiral" # see all possibility in EuclideanExamples.R
samplesize = 100 # number of points on the manifold
SNR = 25 # signal to noise ratio
reg_sampling=TRUE # regular sampling of the point on the manifold or uniformly random
plotTrue= TRUE 

# Generate data
data<- sim_Euclidean_data(name,samplesize,SNR,reg_sampling,plotTrue)

# Estimation of geodesic distances with different methods
Estim<- Geo_estimation(data$true_data,data$discrete_data,data$true_geo,plotTrue,FD_true)

# Assessment of the estimation of a spscific method
mat_to_assess = Estim$estim_geo_true_data # choices are estim_geo_true_data, estim_geo_noisy_data, estim_geo_penalized_isomap, estim_geo_scms 
Rel_err <- assess_goodness_estimation(mat_to_assess,data$true_geo)

#### Analysis of Functional data

FD_true =TRUE

# set up parameters
sce =1 # 1 or 2 (for more info see sim_functional_data.R),K,a,b,SNR,reg_sampling,com_grid,plot_true
K = 30 # number of grid points (each curve is observed on K points on [a,b])
a=-2
b=3
samplesize = 100 # number of points on the manifold
SNR = 0.5 # signal to noise ratio (in Chen and Muller is 0.1 or 0.5)
reg_sampling=TRUE # regular sampling of the point on the manifold or uniformly random
plotTrue= TRUE 
com_grid = 1 # 1 or 0 to indicate if yes or no each curve is observed on a common grid
s = 3 # reduced dimension use for mds and random projection
nb_proj = 100 # number of random projection

# Generate data
data<- sim_functional_data(sce,samplesize,K,a,b,SNR,reg_sampling,com_grid,plotTrue)

# Estimation of geodesic distances with different methods
Estim<- Geo_estimation(data$true_data,data$discrete_data,data$true_geo,plotTrue,FD_true,s,nb_proj,data$grid,data$reg_grid,com_grid)

# Assessment of the estimation of a spscific method
mat_to_assess = Estim$estim_geo_true_data # choices are estim_geo_true_data, estim_geo_noisy_data, estim_geo_smooth_data, estim_geo_penalized_isomap, estim_geo_mds_scms, estim_geo_RP_scms 
Rel_err <- assess_goodness_estimation(mat_to_assess,data$true_geo)


