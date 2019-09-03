
library(reticulate)
library(fda)
library(igraph)
library(caret)
library(vows)
library(KernSmooth)
library(truncnorm)

# this may need to be set at the beginning of the session 
# Susan's
# use_python('/usr/bin/python',required=TRUE)
# use_python('/usr/local/bin/python3.7',required=TRUE)
# use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)

# Marie's
# use_python('/anaconda3/envs/r-reticulate/bin/python',required=TRUE)
use_python('/Users/UQAM/anaconda3/bin/python',required=TRUE)

source('NP_FDA_Ferraty.R')
source('robust_isomap.R')
source('pairwise_geo_estimation.R')
source('pairwise_geo_OUR_tuned.R')

attach(growth)
true_class<- c(rep(1,39),rep(2,54))
true_group<- c(rep(TRUE,39),rep(FALSE,54))
grid<- seq(1,18,length.out=101)
raw_data <- cbind(hgtm,hgtf)
smooth_growth<- smooth_FD_bspline(t(raw_data),age,grid)
detach(growth)


pdf("./Velocity_curves_data_illustration.pdf",width=5,height=8)
par(mar=c(3,4,1,1),mfrow=c(2,1))
matplot(grid,t(smooth_growth[1:39,]),type='l',col="black",lty=1,ylab='Velocity in cm per year')
par(mar=c(4,4,0,1))
matplot(grid,t(smooth_growth[40:93,]),type='l',col="black",lty=2,xlab='Age',ylab='Velocity in cm per year')
dev.off()


# meth <- list("NN" = FALSE,"RD_o" = FALSE,"RD" = FALSE,"SS_o" = FALSE,"SS" = TRUE,"pI" = FALSE,"OUR" = FALSE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE , "L2"=TRUE, "w_L2"= FALSE)
# 
# # estimation of geo distance with methods TRUE in meth
# geo_estim_growth = pairwise_geo_estimation(method=meth,
#                                               true_data=NA,
#                                               discrete_data=smooth_growth,
#                                               true_geo=NA,
#                                               plot_true=FALSE,
#                                               nb_proj=NA,
#                                               grid=NA,
#                                               reg_grid=grid,
#                                               K_dense =length(grid),
#                                               common_grid_true=TRUE,
#                                               Analytic_geo_available=FALSE,
#                                               is_data_smoothed=TRUE)
# 
# 
# saveRDS(geo_estim_growth,file="./geo_estim_vel_growth.rds")
geo_estim_growth<- readRDS("./geo_estim_vel_growth.rds")

# Create balanced train-test partition taking into account the proportion of binary labels in the data
# train.rows[,i] contains the training indices i-th split

# set.seed(73443212)
# nr_train_test_splits = 200
# n_train<- 50
# percentage_train<- n_train/93
# train.rows = createDataPartition(true_group, times = nr_train_test_splits, p = percentage_train, list = FALSE)
# Class_err_discriminant_NP<- list()
# for(geo in 1:length(geo_estim_growth)){
#   Class_err_tmp<- apply(train.rows,2,funopadi.knn.lcv,data=smooth_growth,kind.of.kernel = "quadratic",pairwise_dist=geo_estim_growth[[geo]],true_class=true_class)
#   Class_err_discriminant_NP[[names(geo_estim_growth)[geo]]] <- Class_err_tmp
# }
# write.table(Class_err_discriminant_NP,"./Erreur_classification_velocity.txt",col.names = FALSE,row.names = FALSE)
Class_err<- read.table("./Erreur_classification_velocity.txt")

pdf("./Erreur_classification_velocity.pdf",width=16/4,height=3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.5,1.5,1.5,1.5))
boxplot(Class_err,las=0,names=c('Spline IsoGeo','s=1','s=2','s=3','Spline L2'),main='Classification Error',cex.main=0.75,cex.names=0.75,cex.axis=0.75,cex.lab=0.75)
dev.off()


# ## If we want to compare to other metrics
# for(q_plsr in 5:9){
#   Class_err_tmp<-apply(train.rows,2,funopadi.knn.lcv_ORIGINAL_DISTANCE,data=smooth_growth,kind.of.kernel = "quadratic",semimetric="mplsr",q=q_plsr,true_class=true_class)
#   Class_err_discriminant_NP[[paste("semimetric_mplsr_",q_plsr,sep="")]] <- Class_err_tmp
# }
# 
# for(q_PCA in 4:8){
#   Class_err_tmp<-apply(train.rows,2,funopadi.knn.lcv_ORIGINAL_DISTANCE,data=smooth_growth,kind.of.kernel = "quadratic",semimetric="pca",q=q_PCA,true_class=true_class)
#   Class_err_discriminant_NP[[paste("semimetric_pca_",q_PCA,sep="")]] <- Class_err_tmp
# }
# 
# Class_err_tmp<-apply(train.rows,2,funopadi.knn.lcv_ORIGINAL_DISTANCE,data=smooth_growth,kind.of.kernel = "quadratic",semimetric="L2",true_class=true_class)
# Class_err_discriminant_NP[["L2"]] <- Class_err_tmp
# 
# boxplot(Class_err_discriminant_NP,las=2)


# ########### We can tune h in scms but it is long and give the almost exactly the same results
# Class_err_discriminant_tuned<- list()
# Err_tmp<-rep(0,nr_train_test_splits)
# for(tr in 1:nr_train_test_splits){
#    Err_tmp[tr]<- pairwise_geo_OUR_tuned(smooth_growth,grid,train.rows[,tr],true_class,5)
#  }
#  for(geo in 1:4){
#    Class_err_discriminant_tuned[[paste("tuned_s",geo,sep='')]] <- Err_tmp[,geo]
#  }
# saveRDS(Class_err_discriminant_tuned,file="./Error_classifcation_our_tuned_velocity_growth.rds")

# Class_err_discriminant_tuned<- readRDS("./Error_classifcation_our_tuned_velocity_growth.rds")
# boxplot(Class_err_discriminant_tuned,las=2)

########### knn classification

# kmax=20
# ksNN<- seq(1,kmax,by=2)
# lksNN<-length(ksNN)
# Class_err_discriminant_NP<- list()
# for(geo in 1:length(geo_estim_growth)){
#   Class_err_tmp<- apply(train.rows,2,knn_dist_matrix,group=true_class,ksNN=ksNN,lksNN=lksNN, dist_matrix=geo_estim_growth[[geo]],N=93)
#   Class_err_discriminant_NP[[names(geo_estim_growth)[geo]]] <- Class_err_tmp
# }
# 
# boxplot(Class_err_discriminant_NP,las=2)


########### Clustering based on geo_estim_growth
# 
# Clustering<-fastkmed( geo_estim_growth$estim_geo_OUR3_s2, ncluster = 2, iterate = 100)
# err_tmp<- length(which(Clustering$cluster==true_class))/nrow(smooth_growth)
# min(err_tmp,1-err_tmp)
# 
# Clustering<-fastkmed( geo_estim_growth$estim_geo_OUR3_s3, ncluster = 2, iterate = 100)
# err_tmp<- length(which(Clustering$cluster==true_class))/nrow(smooth_growth)
# min(err_tmp,1-err_tmp)
# 
# Clustering<-fastkmed( geo_estim_growth$estim_geo_OUR3_s4, ncluster = 2, iterate = 100)
# err_tmp<- length(which(Clustering$cluster==true_class))/nrow(smooth_growth)
# min(err_tmp,1-err_tmp)
# 
# Clustering<-fastkmed( geo_estim_growth$estim_geo_SS, ncluster = 2, iterate = 100)
# err_tmp<- length(which(Clustering$cluster==true_class))/nrow(smooth_growth)
# min(err_tmp,1-err_tmp)
# 
# Clustering<-fastkmed( geo_estim_growth$estim_L2, ncluster = 2, iterate = 100)
# err_tmp<- length(which(Clustering$cluster==true_class))/nrow(smooth_growth)
# min(err_tmp,1-err_tmp)



smooth_FD_bspline <- function(discrete_data,grid,grid_dense){
  
  K<- ncol(discrete_data)
  norder<- 6
  knots<- grid
  nbasis <- length(knots) + norder - 2
  #nbasis <- K+2
  basis_obj <- create.bspline.basis(c(grid[1],grid[K]),nbasis, norder, knots)
  loglambda <- seq(-10,-1,by=0.5)
  gcv<-rep(0,length(loglambda))
  for(i in 1:length(loglambda)){
    fd_par_obj <- fdPar(fdobj=basis_obj,Lfdobj=4,lambda=10^loglambda[i])
    smooth_res <- smooth.basis(grid,t(discrete_data),fd_par_obj)
    gcv[i]=mean(smooth_res$gcv)
  }
  ind_m=which(gcv==min(gcv))
  fd_par_obj <- fdPar(fdobj=basis_obj,Lfdobj=4,lambda=10^loglambda[ind_m])

  print(10^loglambda[ind_m])
  smooth_res <- smooth.basis(grid,t(discrete_data),fd_par_obj)
  smooth_curve <- eval.fd(grid_dense,smooth_res$fd,1)
  
  
  
  return(t(smooth_curve))
  
}