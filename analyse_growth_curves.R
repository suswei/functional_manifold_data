#### Analysis of growth data 
## We estimate the geodesic distance matrix with 4 methods
## a) using raw data
## b) using smoothed data
## c) using scms on reduced dim. of smoothed data obtain with random projections
## d) using scms reduced dim. of smoothed data obtain with mds

#### Oracle for number of neigbors and bandwidth h are obtained by minimising the classification errors 

setwd("/Users/marie_uqam/Dropbox/Marie-Moi/Susan_project/manifold learning/Geo_dist_calculation")

library(reticulate)
library(fda)
library(scatterplot3d)
library(kmed)
library(fields)
source('Modif_Isomap.R')
source('Isomap.R')
source('smooth_FD_bspline.R')

# use right version of python
use_python('/anaconda3/bin/python',required=TRUE)
pyIso = import_from_path("getIsomapGdist",path='.')
py_min_neigh = import_from_path("get_min_num_neighbors",path='.')
scms = import_from_path("scms",path='.')

discrete_data <- t(cbind(growth$hgtm,growth$hgtf))
grid <- growth$age
samplesize <- nrow(discrete_data)
K <- length(grid)
true_cluster<- c(rep(1,39),rep(2,54))
nb_proj<-5 # nb projection for RP
s=3 # reduced dimension for mds and random projection


###  a) Direct estimation from the noisy data
num_neigh_min=py_min_neigh$get_min_num_neighbors(discrete_data)
num_neigh=seq(num_neigh_min,samplesize/2,by=2)
misclass_noisy_mani_K<-rep(0,length(num_neigh))
for(j in 1:length(num_neigh)){
  IsomapGdist = pyIso$getIsomapGdist(discrete_data,num_neigh[j])
  res<-fastkmed( IsomapGdist, ncluster = 2, iterate = 50)
  temp=table(res$cluster,true_cluster)
  misclass_noisy_mani_K[j]=min(1-(sum(diag(temp))/samplesize),sum(diag(temp))/samplesize)
}
err_noisy=min(misclass_noisy_mani_K)
ind_op=min(which(misclass_noisy_mani_K==err_noisy))
Isomap_R = Isomap(discrete_data,2,num_neigh[ind_op])
Iso_data <- Isomap_R[[1]][[1]]
res<-fastkmed( Isomap_R[[2]], ncluster = 2, iterate = 50)

par(mfrow=c(2,2))
matplot(grid,t(discrete_data),main="Observed data",type='l', col=c(rep('blue',39),rep('red',54))) 
image.plot(Isomap_R[[2]],main='geo raw data')
plot(Iso_data[,1],Iso_data[,2],type='n',main="Mani comp.")
points(Iso_data[,1],Iso_data[,2],col=c(rep('blue',39),rep('red',54)))
plot(Iso_data[,1],Iso_data[,2],type='n',main=paste('misclass. raw data = ',signif(err_noisy,digits=4),sep=''))
points(Iso_data[,1],Iso_data[,2],col=res$cluster)


### b) Smooth the data with bspline + roughness penalty and direct estimation
smooth_data = smooth_FD_bspline(discrete_data,matrix(grid,ncol=K,nrow=samplesize,byrow=TRUE),grid)

num_neigh_min=py_min_neigh$get_min_num_neighbors(smooth_data)
num_neigh=seq(num_neigh_min,samplesize/2,by=2)
misclass_smooth_mani_K= rep(0,length(num_neigh))
for(j in 1:length(num_neigh)){
  IsomapGdist = pyIso$getIsomapGdist(smooth_data,num_neigh[j])
  res<-fastkmed( IsomapGdist, ncluster = 2, iterate = 50)
  temp=table(res$cluster,true_cluster)
  misclass_smooth_mani_K[j]=min(1-(sum(diag(temp))/samplesize),sum(diag(temp))/samplesize)
}
err_smooth=min(misclass_smooth_mani_K)
ind_op=min(which(misclass_smooth_mani_K==err_smooth))
Isomap_R = Isomap(smooth_data,2,num_neigh[ind_op])
Iso_data <- Isomap_R[[1]][[1]]
res<-fastkmed( Isomap_R[[2]], ncluster = 2, iterate = 50)

par(mfrow=c(2,2))
matplot(grid,t(smooth_data),main="Smooth data",type='l', col=c(rep('blue',39),rep('red',54))) 
image.plot(Isomap_R[[2]],main='geo smooth data')
plot(Iso_data[,1],Iso_data[,2],type='n',main="Mani comp.")
points(Iso_data[,1],Iso_data[,2],col=c(rep('blue',39),rep('red',54)))
plot(Iso_data[,1],Iso_data[,2],type='n',main=paste('misclass. raw data = ',signif(err_smooth,digits=4),sep=''))
points(Iso_data[,1],Iso_data[,2],col=res$cluster)


### c) Smooth the data with bspline + roughness penalty use random projection to obtain a s-dimensional 
# version of each curve and use SCMS to estimate the manifold.

scms_h=seq(0.1,2,by=0.2)
op_h_proj<-rep(0,nb_proj)
op_num_nei_proj<-rep(0,nb_proj)
err_proj<-rep(0,nb_proj)
geo_proj<-matrix(0,nrow=nb_proj,ncol=samplesize*(samplesize-1)/2)
for(proj in 1:nb_proj){
  randprojmat = matrix(rnorm(s*K,mean=0,sd=1/sqrt(s)), nrow=s, ncol=K) 
  projected = randprojmat %*% t(smooth_data)
  Error_denoised_h= matrix(0,nrow=length(scms_h),ncol=2)
  for(bw in 1:length(scms_h)){
    denoised = scms$scms(t(projected), scms_h[bw])
    num_neigh_min=py_min_neigh$get_min_num_neighbors(denoised)
    num_neigh=seq(num_neigh_min,samplesize/2,by=4)
    Error_denoised_K= rep(0,length(num_neigh))
    for(j in 1:length(num_neigh)){
      IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,num_neigh[j])
      res<-fastkmed( IsomapGdist_denoised, ncluster = 2, iterate = 50)
      temp=table(res$cluster,true_cluster)
      Error_denoised_K[j]=min(1-(sum(diag(temp))/samplesize),sum(diag(temp))/samplesize)
    }
    ind_neigh=min(which(Error_denoised_K==min(Error_denoised_K)))
    Error_denoised_h[bw,]=c(Error_denoised_K[ind_neigh],num_neigh[ind_neigh])
  }
  ind=min(which(Error_denoised_h[,1]==min(Error_denoised_h[,1])))
  denoised = scms$scms(t(projected), scms_h[ind])
  IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,Error_denoised_h[ind,2])
  geo_proj[proj,]<-IsomapGdist_denoised[lower.tri(IsomapGdist_denoised, diag = FALSE)]
  op_h_proj[proj]=scms_h[ind]
  op_num_nei_proj[proj]=Error_denoised_h[ind,2]
  err_proj[proj]=min(Error_denoised_h[,1])
}
 
ensemble_geo<-apply(geo_proj,2,mean,trim=0.25)
ens_geo_mat<-matrix(0,samplesize,samplesize)
ens_geo_mat[lower.tri(ens_geo_mat, diag = FALSE)]=ensemble_geo
ens_geo_mat= ens_geo_mat + t(ens_geo_mat)

res<-fastkmed( ens_geo_mat, ncluster = 2, iterate = 50)
temp=table(res$cluster,true_cluster)
err=min(1-(sum(diag(temp))/samplesize),sum(diag(temp))/samplesize)
Iso_data=Modif_Isomap(ens_geo_mat,2)

par(mfrow=c(1,3))
image.plot(Isomap_R[[2]],main='geo RP + scms ensemble')
plot(Iso_data[[1]][,1],Iso_data[[1]][,2],xlab='',ylab='',type='n',main='Mani comp')
points(Iso_data[[1]][,1],Iso_data[[1]][,2],col=c(rep('blue',39),rep('red',54)))
plot(Iso_data[[1]][,1],Iso_data[[1]][,2],xlab='',ylab='',type='n',main=paste('misclass. RP = ',signif(err,digits=4),sep=''))
points(Iso_data[[1]][,1],Iso_data[[1]][,2],col=res$cluster)


### d) Smooth the data with bspline + roughness penalty use random projection to obtain a s-dimensional 
# version of each curve and use SCMS to estimate the manifold.

euc_dis_mat=dist(discrete_data)
mds_mat = cmdscale(euc_dis_mat,eig=TRUE, k=s) 
projected = mds_mat$points

scms_h=seq(0.1,2,by=0.2)
Error_denoised_h= matrix(0,nrow=length(scms_h),ncol=2)
for(bw in 1:length(scms_h)){
  denoised = scms$scms(projected, scms_h[bw])
  num_neigh_min=py_min_neigh$get_min_num_neighbors(denoised)
  num_neigh=seq(num_neigh_min,samplesize/2,by=4)
  Error_denoised_K= rep(0,length(num_neigh))
  for(j in 1:length(num_neigh)){
    IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,num_neigh[j])
    res<-fastkmed( IsomapGdist_denoised, ncluster = 2, iterate = 50)
    temp=table(res$cluster,true_cluster)
    Error_denoised_K[j]=min(1-(sum(diag(temp))/samplesize),sum(diag(temp))/samplesize)
  }
  ind_neigh=min(which(Error_denoised_K==min(Error_denoised_K)))
  Error_denoised_h[bw,]=c(Error_denoised_K[ind_neigh],num_neigh[ind_neigh])
}
err_mds = min(Error_denoised_h[,1])
ind=min(which(Error_denoised_h[,1]==err_mds))
denoised = scms$scms(projected, scms_h[ind])
Isomap_R = Isomap(denoised,2,Error_denoised_h[ind,2])
Iso_data <- Isomap_R[[1]][[1]]
res<-fastkmed( Isomap_R[[2]], ncluster = 2, iterate = 50)

par(mfrow=c(1,3))
image.plot(Isomap_R[[2]],main='geo mds + scms data')
plot(Iso_data[,1],Iso_data[,2],type='n',main="Mani comp.")
points(Iso_data[,1],Iso_data[,2],col=c(rep('blue',39),rep('red',54)))
plot(Iso_data[,1],Iso_data[,2],type='n',main=paste('misclass. mds = ',signif(err_mds,digits=4),sep=''))
points(Iso_data[,1],Iso_data[,2],col=res$cluster)


