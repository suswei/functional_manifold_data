# Title     : compare_methods.R
# Objective : This is modeled after Marie's main.R but designed for submission to Spartan HPC. Performs sweep of certain parameters.
# Created by: suswei
# Created on: 3/12/18

# library(DescTools)
library(reticulate)
library(fields)
library(fda)
library(matlabr)
library(igraph)

source('EuclideanExamples.R')
source('full_geo_from_adj_geo.R')
source('sim_functional_data.R')
source('Geo_estimation.R')
source('pairwise_geo_estimation.R')
source('sim_Euclidean_data.R')
source('assess_goodness_estimation.R')
source('robust_isomap.R')

# get sweeping id
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
slurm_arrayid = as.numeric(slurm_arrayid)
print(slurm_arrayid)

# parameters that will be fixed during the sim study
K = 30 # number of grid points (each curve is observed on K points on [a,b])
com_grid = 1 # 1 or 0 to indicate if yes or no each curve is observed on a common grid
plotTrue = FALSE
FD_true = TRUE

# parameters under study
# scenario
# samplesize # number of points on the manifold
# SNR # signal to noise ratio (in Chen and Muller is 0.1 or 0.5)
# reg_sampling {True,False} # regular sampling of the point on the manifold or uniformly random

sces = c(1,2,4)
samplesizes = c(100,250)
SNRs = c(0.1,0.5)
reg_samplings = c(TRUE,FALSE)
mcs = 1:100

unravel=arrayInd(slurm_arrayid,c(length(sces),length(samplesizes), length(SNRs), length(reg_samplings), length(mcs)))

# actual parameters for this run
sce = sces[unravel[1,1]]
samplesize = samplesizes[unravel[1,2]]
SNR = SNRs[unravel[1,3]]
reg_sampling = reg_samplings[unravel[1,4]]
mc = mcs[unravel[1,5]]


# Generate data
data<- sim_functional_data(sce, samplesize, K, SNR, reg_sampling, com_grid, plotTrue)

# Estimation of geodesic distances with different methods
meth <- list("NN" = TRUE,"RD_o" = TRUE,"RD" = TRUE,"SS_o" = TRUE,"SS" = TRUE,"pI" = FALSE,"OUR" = TRUE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE )# see pairwise_geo_estimation for more info
Estim<- pairwise_geo_estimation(meth,data$noiseless_data,data$noisy_data,data$analytic_geo,plotTrue,FD_true,nb_proj,data$grid,data$reg_grid,com_grid)
# TODO: decide if we need to save Estim?
# saveRDS(Estim, file = sprintf("sce=%d_samplesize=%d_SNR=%d_reg_sampling=%d_s=%d_mc=%d",sce,samplesize,SNR,reg_sampling,s,mc))

# Assessment of the estimation of a spscific method
Rel_errs = lapply(Estim, assess_goodness_estimation, true_geo = data$analytic_geo)
print(Rel_errs)

saveRDS(Rel_errs, file = sprintf("taskid=%d_sce=%d_samplesize=%d_SNR=%.01f_reg_sampling=%d_mc=%d",slurm_arrayid,sce,samplesize,SNR,reg_sampling,mc))
# saveRDS(Rel_errs, file = sprintf("/data/cephfs/punim0715/taskid=%d_sce=%d_samplesize=%d_SNR=%.01f_reg_sampling=%d_s=%d_mc=%d",slurm_arrayid,sce,samplesize,SNR,reg_sampling,s,mc))

# open with readRDS(filename)