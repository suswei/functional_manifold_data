# Title     : TODO
# Objective : This is modeled after Marie's main.R. Simulation sweep of parameters.
# Created by: suswei
# Created on: 3/12/18

# parameters that will be fixed during the sim study
# K = 30 # number of grid points (each curve is observed on K points on [a,b])
# com_grid = 1 # 1 or 0 to indicate if yes or no each curve is observed on a common grid
# nb_proj = 500 # number of random projection
# plotTrue = FALSE
# a = -2
# b = 3

# parameters under study
# sce {1,2} (for more info see sim_functional_data.R),K,a,b,SNR,reg_sampling,com_grid,plot_true
# samplesize # number of points on the manifold
# SNR # signal to noise ratio (in Chen and Muller is 0.1 or 0.5)
# reg_sampling=TRUE # regular sampling of the point on the manifold or uniformly random
# s = 3 # reduced dimension use for mds and random projection

# spartanSim <- function(sce, samplesize, SNR, reg_sampling, s, mc, K = 30, com_grid=1, nb_proj = 500, plotTrue = FALSE, a = -2, b = 3, FD_true = TRUE){
spartanSim <- function(mc, K = 30, com_grid=1, nb_proj = 500, plotTrue = FALSE, a = -2, b = 3, FD_true = TRUE){

sce = 1
samplesize = 100
SNR = 0.5
reg_sampling = TRUE
s = 2

rm(list = ls())

source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('Modif_Isomap.R')
source('pairwiseDistances.R')
source('full_geo_from_adj_geo.R')
source('smooth_FD_bspline.R')
source('sim_functional_data.R')
source('Geo_estimation.R')
source('sim_Euclidean_data.R')
source('assess_goodness_estimation.R')

library(fields)
library(reticulate)
library(fda)
library(matlabr)

pyIso = import_from_path("getIsomapGdist",path='.')
py_min_neigh = import_from_path("get_min_num_neighbors",path='.')
scms = import_from_path("scms",path='.')

# Generate data
data<- sim_functional_data(sce,samplesize,K,a,b,SNR,reg_sampling,com_grid,plotTrue)

# Estimation of geodesic distances with different methods
Estim<- Geo_estimation(data$true_data,data$discrete_data,data$true_geo,plotTrue,FD_true,s,nb_proj,data$grid,data$reg_grid,com_grid)

saveRDS(Estim, file = paste("sce=%d,samplesize=%d,SNR=%d,reg_sampling=%d,s=%d,mc=%d",sce,samplesize,SNR,reg_sampling,s,mc))

# Assessment of the estimation of a spscific method
# mat_to_assess = Estim$estim_geo_true_data # choices are estim_geo_true_data, estim_geo_noisy_data, estim_geo_smooth_data, estim_geo_penalized_isomap, estim_geo_mds_scms, estim_geo_RP_scms
# Rel_err <- assess_goodness_estimation(mat_to_assess,data$true_geo)

}