rm(list = ls())

setwd("/Users/marie_uqam/Dropbox/Marie-Moi/Susan_project/manifold learning/Geo_dist_calculation")


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
scms = import_from_path("scms",path='.')

# set up parameters

#name = c("archimedean-spiral",'circle',"right-angle","sin-cos-curve")
name = c('circle',"sin-cos-curve")
samplesize = c(200,500)
#samplesize = c(500)
SNR = c(20,25,30) # noise sd
#SNR = c(20) # noise sd
#reg_sampling=c(TRUE,FALSE)
reg_sampling=c(FALSE)



#Obtain the plot for each combination of parameter

for(mani_type in 1:length(name)){
  for(n in 1:length(samplesize)){
    for(level_noise in 1:length(SNR)){
      for(sampling in 1:length(reg_sampling)){
        sim_study(name[mani_type],samplesize[n],SNR[level_noise],reg_sampling[sampling])
      }
    }
  }
}
