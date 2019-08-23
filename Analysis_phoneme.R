# Classification phoneme data, comparison with results from Delaigle and Hall

library(reticulate)
library(fda)
library(scatterplot3d)
library(igraph)
library(ROCR)
library(caret)
library(vows)
library(kmed)
library(ElemStatLearn)
library(KernSmooth)

# this may need to be set at the beginning of the session 
# Susan's
# use_python('/usr/bin/python',required=TRUE)
# use_python('/usr/local/bin/python3.7',required=TRUE)
# use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)

# Marie's
# use_python('/anaconda3/envs/r-reticulate/bin/python',required=TRUE)
use_python('/Users/UQAM/anaconda3/bin/python',required=TRUE)


source('robust_isomap.R')
source('pairwise_geo_estimation.R')
source("weight_L2.R")
source("processRealData.R")
source("NP_FDA_Ferraty.R")

PhonData_aa<-phoneme[which(phoneme$g=="aa"),]
PhonData_ao<-phoneme[which(phoneme$g=="ao"),]
Phon_aa_ao_temp<- rbind(PhonData_aa,PhonData_ao)
n_aa<-nrow(PhonData_aa)
n_ao<-nrow(PhonData_ao)
n<- n_aa+n_ao
true_group<- rep(FALSE,n)
true_group[which(Phon_aa_ao_temp$g=="aa")] <-TRUE
true_class<- rep(0,n)
true_class[true_group]<- 1
Phon_aa_ao<- as.matrix(Phon_aa_ao_temp[,-c(257,258)])
grid<- 1:256

# We smooth the data with local linear smoother
Smooth_Phon<-matrix(nrow=n,ncol=256)
for(i in 1:n){
  Smooth_Phon[i,]<-locpoly(grid,Phon_aa_ao[i,],degree=1,kernel="normal",bandwidth=4,gridsize=256)$y
}
matplot(grid,t(Smooth_Phon[1:n_aa,]),type='l')
matplot(grid,t(Smooth_Phon[(n_aa+1):n,]),type='l')

# we just considerd grid = [1,50] as in Delaigle and Hall and we take a subset of each group (dataset too big)
smooth_phon_trunc<-  Smooth_Phon[c(1:250,(n_aa+1):(n_aa+250)),1:50]
true_group_trunc<- true_group[c(1:250,(n_aa+1):(n_aa+250))]
true_class_trunc<- true_class[c(1:250,(n_aa+1):(n_aa+250))]
grid_trunc<- 1:50


meth <- list("NN" = FALSE,"RD_o" = FALSE,"RD" = FALSE,"SS_o" = FALSE,"SS" = TRUE,"pI" = FALSE,"OUR" = FALSE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE , "L2"=TRUE, "w_L2"= FALSE)

geo_estim = pairwise_geo_estimation(method=meth,
                                    true_data=NA,
                                    discrete_data=smooth_phon_trunc,
                                    true_geo=NA,
                                    plot_true=FALSE,
                                    nb_proj=NA,
                                    grid=matrix(rep(grid,length(grid_trunc)),nrow=length(grid_trunc),byrow = TRUE),
                                    reg_grid=grid_trunc,
                                    K_dense =100,
                                    common_grid_true=TRUE,
                                    Analytic_geo_available=FALSE,
                                    is_data_smoothed=TRUE)
write.table(geo_estim,file="./geo_estim_phoneme.txt",col.names = FALSE, row.names = FALSE)
# Assess classificaton performance
nr_train_test_splits = 200
n_train<- 100
percentage_train<- n_train/nrow(smooth_phon_trunc)
# Create balanced train-test partition taking into account the proportion of binary labels in the data
# train.rows[,i] contains the training indices i-th split
train.rows = createDataPartition(true_group_trunc, times = nr_train_test_splits, p = percentage_train, list = FALSE)

Class_err_discriminant_NP<- list()
for(geo in 1:length(geo_estim)){
  Class_err_tmp<- apply(train.rows,2,funopadi.knn.lcv,data=smooth_phon_trunc,kind.of.kernel = "quadratic",pairwise_dist=geo_estim[[geo]],true_class=true_class_trunc+1)
  Class_err_discriminant_NP[[names(geo_estim)[geo]]] <- Class_err_tmp
}
Class_err_discriminant_NP_ori<- apply(train.rows,2,funopadi.knn.lcv_ORIGINAL_DISTANCE,data=smooth_phon_trunc,true_class=true_class_trunc+1,kind.of.kernel = "quadratic",q=2,nknot=10,range.grid=c(1,50))
Class_err_discriminant_NP[["semimetric_Ferraty"]]<-Class_err_discriminant_NP_ori
boxplot(Class_err_discriminant_NP,las=2,ylab="misclassification rate",names=c(substr(names(geo_estim)[1:4], 11, 17),substr(names(geo_estim)[5], 7, 10),"semimetric_NP"),main="Classification error funopadi.knn.lcv")

Class_err_NP_reg<- list()
for(geo in 1:length(geo_estim)){
  Class_err_tmp<- apply(train.rows,2,funopare.knn.gcv,data=smooth_phon_trunc,kind.of.kernel = "quadratic",pairwise_dist=geo_estim[[geo]],true_class=true_class_trunc)
  Class_err_NP_reg[[names(geo_estim)[geo]]] <- Class_err_tmp
}
Class_err__NP_reg_ori<- apply(train.rows,2,funopare.knn.gcv_ORIGINAL_DISTANCE,data=smooth_phon_trunc,true_class=true_class_trunc,kind.of.kernel = "quadratic",q=2,nknot=10,range.grid=c(1,50))
Class_err_NP_reg[["semimetric_Ferraty"]]<-Class_err__NP_reg_ori

boxplot(Class_err_NP_reg,las=2,ylab="misclassification rate",names=c(substr(names(geo_estim)[1:4], 11, 17),substr(names(geo_estim)[5], 7, 10),"semimetric_NP"),main="Classification error funopare.knn.gcv")

