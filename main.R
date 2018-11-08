rm(list = ls())

source('EuclideanExamples.R')
source('functions.R')
source('Isomap.R')
source('pairwiseDistances.R')
source('full_geo_from_adj_geo.R')
source('sim_study.R')

library(fields)
library(reticulate)


# use right version of python
use_python('/anaconda3/bin/python',required=TRUE)
pyIso = import_from_path("getIsomapGdist",path='.')
pyIsoAuto = import_from_path("get_auto_isomap_gdist",path='.')

# set up parameters
name = c("archimedean-spiral",'circle',"right-angle","sin-cos-curve")
samplesize = c(100,200,500)
noi_sd = c(0.1,0.5,1,1.5) # noise sd
reg_sampling=c(TRUE,FALSE)


#Obtain the plot for each combination of parameter

for(mani_type in 1:length(name)){
  for(n in 1:length(samplesize)){
    for(level_noise in 1:length(noi_sd)){
      for(sampling in 1:2){
        sim_study(name[mani_type],samplesize[n],noi_sd[level_noise],reg_sampling[sampling])
      }
    }
  }
}
