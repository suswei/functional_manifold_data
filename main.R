library(reticulate)
# use right version of python
# use_python('/Users/UQAM/anaconda3/bin/python',required=TRUE)
use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)

rm(list = ls())

source('EuclideanExamples.R')
source('full_geo_from_adj_geo.R')
source('sim_functional_data.R')
source('Geo_estimation.R')
source('pairwise_geo_estimation.R')
source('sim_Euclidean_data.R')
source('assess_goodness_estimation.R')

library(DescTools)
library(fields)
library(fda)
library(matlabr)


pyIso = import_from_path("getIsomapGdist",path='.')
py_min_neigh = import_from_path("get_min_num_neighbors",path='.')
scms = import_from_path("scms",path='.')


#### Analysis of Functional data

FD_true =TRUE

# set up parameters
sce =3 # 1, 2 or 3 (for more info see sim_functional_data.R)
K = 30 # number of grid points (each curve is observed on K points on [a,b])
samplesize = 100 # number of points on the manifold
SNR = 0.1 # signal to noise ratio (in Chen and Muller is 0.1 or 0.5)
reg_sampling=TRUE # regular sampling of the point on the manifold or uniformly random
plotTrue= TRUE 
com_grid = 1 # 1 or 0 to indicate if yes or no each curve is observed on a common grid
nb_proj = 20 # number of random projection
meth<- list("NN" = TRUE,"RD" = TRUE,"SS" = TRUE,"pI" = FALSE,"OUR" = TRUE,"RP" = TRUE ) # see pairwise_geo_estimation for more info

# Generate data
data<- sim_functional_data(sce,samplesize,K,SNR,reg_sampling,com_grid,plotTrue)

# Estimation of geodesic distances with different methods
Estim<- pairwise_geo_estimation(meth,data$noiseless_data,data$noisy_data,data$analytic_geo,plotTrue,FD_true,nb_proj,data$grid,data$reg_grid,com_grid)

# Comparision of the different methods
nom_methode <- rep(names(Estim),rep(3,length(Estim)))
measure <- rep(c("rmse","epsilon_isometry_auc","pearson_corr"),length(Estim))
resu<- rep(0,3*length(Estim))
res<-data.frame(nom_methode,measure,resu)
res$measure<- as.factor(res$measure)
for(i in 1:length(Estim)){
  res[(1+(i-1)*3):(3*i),3]<- unlist(assess_goodness_estimation(Estim[[i]],data$analytic_geo))
}
par(mfrow=c(1,1))
interaction.plot(res$measure,res$nom_methode,res$resu, type = "b", col=1:6,xlab="different measures",ylab="different methods")

# Assessment of the estimation of a specific method
mat_to_assess = Estim$estim_geo_true_data # choices are estim_geo_true_data, estim_geo_noisy_data, estim_geo_smooth_data, estim_geo_penalized_isomap, estim_geo_mds_scms, estim_geo_RP_scms 
Rel_err <- assess_goodness_estimation(mat_to_assess,data$analytic_geo)
print(Rel_err)



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