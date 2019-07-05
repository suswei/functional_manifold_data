## Analysis of Growth curves using the smoothing strategy of Chen and al in Biometrics

library(reticulate)
library(fda)
library(scatterplot3d)
library(igraph)
library(ROCR)
library(caret)
library(vows)

# this may need to be set at the beginning of the session 
# Susan's
# use_python('/usr/bin/python',required=TRUE)
# use_python('/usr/local/bin/python3.7',required=TRUE)
# use_python('/Users/suswei/anaconda3/bin/python',required=TRUE)

# Marie's
# use_python('/anaconda3/envs/r-reticulate/bin/python',required=TRUE)
# use_python('/Users/UQAM/anaconda3/bin/python',required=TRUE)

source("weight_L2.R")
source('robust_isomap.R')
source('pairwise_geo_estimation.R')
source("processRealData.R")
source("NP_FDA_Ferraty.R")


attach(growth)
true_class<- c(rep(1,39),rep(2,54))
true_group<- c(rep(TRUE,39),rep(FALSE,54))
reg_grid<- seq(1,18,length.out=101)

raw_data <- t(cbind(hgtm,hgtf))
temp<- pairwise_weigthed_L2(raw_data,age,reg_grid)
smooth_data_ori<- temp$smooth_data
estim_weigthed_L2_ori <- temp$pairwise_L2

smooth_vel<- (smooth_data_ori[,2:ncol(smooth_data_ori)] - smooth_data_ori[,1:(ncol(smooth_data_ori)-1)])/ (reg_grid[2]-reg_grid[1])
temp2<- pairwise_weigthed_L2(smooth_vel,reg_grid[-1],reg_grid[-1])
smooth_data_vel<- temp2$smooth_data
estim_weigthed_L2_vel <- temp2$pairwise_L2
detach(growth)


#### Classification of velocity curve
# We compare results with SS, OUR3, L2 and weight L2 

meth <- list("NN" = FALSE,"RD_o" = FALSE,"RD" = FALSE,"SS_o" = FALSE,"SS" = TRUE,"pI" = FALSE,"OUR" = FALSE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE , "L2"=TRUE, "w_L2"= FALSE)
Class_err<- list()

# estimation of geo distance with methods TRUE in meth
geo_estim = pairwise_geo_estimation(method=meth,
                                      true_data=NA,
                                      discrete_data=smooth_data_vel,
                                      true_geo=NA,
                                      plot_true=FALSE,
                                      nb_proj=NA,
                                      grid=NA,
                                      reg_grid=reg_grid[-1],
                                      K_dense =100,
                                      common_grid_true=TRUE,
                                      Analytic_geo_available=FALSE,
                                      is_data_smoothed=TRUE)
  
# Assess classificaton performance
nr_train_test_splits = 500
  
# Create balanced train-test partition taking into account the proportion of binary labels in the data
# train.rows[,i] contains the training indices i-th split
train.rows = createDataPartition(true_group, times = nr_train_test_splits, p = 0.53, list = FALSE)
  
## Weigthed L2
Class_err_tmp<- apply(train.rows,2,funopadi.knn.lcv,data=smooth_data_vel,kind.of.kernel = "quadratic",pairwise_dist=estim_weigthed_L2_vel,true_class=true_class)
Class_err[["weigthed_L2"]] <- Class_err_tmp
  
# Other methods
for(geo in 1:length(geo_estim)){
  Class_err_tmp<- apply(train.rows,2,funopadi.knn.lcv,data=smooth_data_vel,kind.of.kernel = "quadratic",pairwise_dist=geo_estim[[geo]],true_class=true_class)
  Class_err[[names(geo_estim)[geo]]] <- Class_err_tmp
}
  
# Boxplot

pdf("classification_results_velocity_curve.pdf",width=12,height=5)
boxplot(Class_err,las=2,ylab="misclassification rate",names=c("Weigthed_L2",substr(names(geo_estim)[1:4], 11, 17),substr(names(geo_estim)[5], 7, 10)),main="Classification error for velocity curves")
dev.off()



#### Classification of growth curve
# We compare results with SS, OUR3, L2 and weight L2 

meth <- list("NN" = FALSE,"RD_o" = FALSE,"RD" = FALSE,"SS_o" = FALSE,"SS" = TRUE,"pI" = FALSE,"OUR" = FALSE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE , "L2"=TRUE, "w_L2"= FALSE)
Class_err<- list()

# estimation of geo distance with methods TRUE in meth
geo_estim = pairwise_geo_estimation(method=meth,
                                    true_data=NA,
                                    discrete_data=smooth_data_ori,
                                    true_geo=NA,
                                    plot_true=FALSE,
                                    nb_proj=NA,
                                    grid=NA,
                                    reg_grid=reg_grid,
                                    K_dense =100,
                                    common_grid_true=TRUE,
                                    Analytic_geo_available=FALSE,
                                    is_data_smoothed=TRUE)

# Assess classificaton performance
nr_train_test_splits = 500

# Create balanced train-test partition taking into account the proportion of binary labels in the data
# train.rows[,i] contains the training indices i-th split
train.rows = createDataPartition(true_group, times = nr_train_test_splits, p = 0.53, list = FALSE)

## Weigthed L2
Class_err_tmp<- apply(train.rows,2,funopadi.knn.lcv,data=smooth_data_vel,kind.of.kernel = "quadratic",pairwise_dist=estim_weigthed_L2_vel,true_class=true_class)
Class_err[["weigthed_L2"]] <- Class_err_tmp

# Other methods
for(geo in 1:length(geo_estim)){
  Class_err_tmp<- apply(train.rows,2,funopadi.knn.lcv,data=smooth_data_vel,kind.of.kernel = "quadratic",pairwise_dist=geo_estim[[geo]],true_class=true_class)
  Class_err[[names(geo_estim)[geo]]] <- Class_err_tmp
}

# Boxplot

pdf("classification_results_growth_curve.pdf",width=12,height=5)
boxplot(Class_err,las=2,ylab="misclassification rate",names=c("Weigthed_L2",substr(names(geo_estim)[1:4], 11, 17),substr(names(geo_estim)[5], 7, 10)),main="Classification error for growth curves")
dev.off()

