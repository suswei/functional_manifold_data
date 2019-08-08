# Title     : compare_methods.R
# Objective : compare multiple estimators of geodesic distance on simulated data. this is called by compare_methods.sh!
# Created by: suswei
# Created on: 3/12/18

# library(DescTools)
library(reticulate)
use_python("/Users/suswei/anaconda3/bin/python",required=TRUE)
library(fields)
library(fda)
library(igraph)
library(matlabr)

source('EuclideanExamples.R')
source('full_geo_from_adj_geo.R')
source('sim_functional_data.R')
source('Geo_estimation.R')
source('pairwise_geo_estimation.R')
source('sim_Euclidean_data.R')
source('assess_goodness_estimation.R')
source('robust_isomap.R')



# parameters that will be fixed during the sim study
com_grid = 1 # 1 or 0 to indicate if yes or no each curve is observed on a common grid
plotTrue = FALSE
samplesize = 100
Ks_smooth = 100

## STUDY THE EFFECT OF THE FOLLOWING PARAMETERS 
mcs = 1:10
# parameters in sim_functional_data
sces = c(5,2,4)
SNRs = c(0.1,0.5)
reg_samplings = c(TRUE,FALSE)
Ks_obs = c(100)

total_tasks = length(sces)*length(SNRs)*length(reg_samplings)*length(Ks_obs)*length(mcs)


for(slurm_arrayid in 1:total_tasks){
  
  # Hardcode
  unravel=arrayInd(slurm_arrayid,c(length(sces), length(SNRs), length(reg_samplings), length(Ks_obs),length(mcs)))
  
  # actual parameters for this run
  sce = sces[unravel[1,1]]
  SNR = SNRs[unravel[1,2]]
  reg_sampling = reg_samplings[unravel[1,3]]
  K_obs = Ks_obs[unravel[1,4]]
  mc = mcs[unravel[1,5]]
  
  
  # Generate data
  data<- sim_functional_data(samplesize = samplesize, 
                             com_grid = com_grid, 
                             plot_true = plotTrue,
                             ### parameters under study
                             sce = sce, 
                             SNR = SNR, 
                             reg_sampling = reg_sampling,
                             K = K_obs)
  
  
  # Estimation of geodesic distances with different methods
  meth <- list("NN" = FALSE,
               "RD_o" = FALSE,
               "RD" = FALSE,
               "SS_o" = FALSE,
               "SS" = FALSE,
               "pI" = TRUE,
               "OUR" = FALSE,
               "OUR2" = FALSE,
               "OUR3"=FALSE,
               "RP" = FALSE,
               "L2" = FALSE, 
               "w_L2" = FALSE )# see pairwise_geo_estimation for more info
  
  Estim<- pairwise_geo_estimation(method=meth,
                                  true_data = data$noiseless_data,
                                  discrete_data = data$noisy_data,data$analytic_geo,
                                  plot_true = plotTrue,
                                  grid = data$grid,
                                  reg_grid = data$reg_grid,
                                  common_grid_true =com_grid,
                                  K_dense = Ks_smooth)
}
