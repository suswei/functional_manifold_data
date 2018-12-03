# Title     : TODO
# Objective : This is modeled after Marie's main.R. Simulation sweep of parameters.
# Created by: suswei
# Created on: 3/12/18

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
slurm_arrayid = as.numeric(slurm_arrayid)
print(slurm_arrayid)

# parameters that will be fixed during the sim study
K = 30 # number of grid points (each curve is observed on K points on [a,b])
com_grid = 1 # 1 or 0 to indicate if yes or no each curve is observed on a common grid
nb_proj = 100 # number of random projection
plotTrue = FALSE
a = -2
b = 3
FD_true = TRUE

# parameters under study
# sce {1,2} (for more info see sim_functional_data.R),K,a,b,SNR,reg_sampling,com_grid,plot_true
# samplesize # number of points on the manifold
# SNR # signal to noise ratio (in Chen and Muller is 0.1 or 0.5)
# reg_sampling=TRUE # regular sampling of the point on the manifold or uniformly random
# s = {1,2,3,4,5} # reduced dimension use for mds and random projection

sces = c(1,2)
samplesizes = c(100,250)
SNRs = c(0.1,0.5)
reg_samplings = c(TRUE,FALSE)
ss = c(1,2,3,4)
mcs = 1:100

unravel=arrayInd(slurm_arrayid,c(length(sces), length(samplesizes), length(SNRs), length(reg_samplings), length(ss), length(mcs)))

# actual parameters for this run
sce = sces[unravel[1,1]]
samplesize = samplesizes[unravel[1,2]]
SNR = SNRs[unravel[1,3]]
reg_sampling = reg_samplings[unravel[1,4]]
s = ss[unravel[1,5]]
mc = mcs[unravel[1,6]]

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

# TODO: it doesn't make ansy sense to define these outside of Geo_estimation?
pyIso = import_from_path("getIsomapGdist",path='.')
py_min_neigh = import_from_path("get_min_num_neighbors",path='.')
scms = import_from_path("scms",path='.')

# Generate data
data<- sim_functional_data(sce,samplesize,K,a,b,SNR,reg_sampling,com_grid,plotTrue)

# Estimation of geodesic distances with different methods
Estim<- Geo_estimation(data$true_data,data$discrete_data,data$true_geo,plotTrue,FD_true,s,nb_proj,data$grid,data$reg_grid,com_grid)

# TODO: decide if we need to save Estim?
# saveRDS(Estim, file = sprintf("sce=%d_samplesize=%d_SNR=%d_reg_sampling=%d_s=%d_mc=%d",sce,samplesize,SNR,reg_sampling,s,mc))

# Assessment of the estimation of a spscific method
Rel_errs = lapply(Estim,assess_goodness_estimation, true_geo = data$true_geo)


saveRDS(Rel_errs, file = sprintf("sce=%d_samplesize=%d_SNR=%.01f_reg_sampling=%d_s=%d_mc=%d",sce,samplesize,SNR,reg_sampling,s,mc))
