#### Classification of real dataset preprocessed by processRealData.R
# We compare results with SS, OUR3, L2 and weight L2 

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


source('robust_isomap.R')
source('pairwise_geo_estimation.R')
source("weight_L2.R")
source("processRealData.R")



# call processRealData which each element is a list containing "discrete_data" , "reg_grid"  "true_group"
dataset = processRealData()

meth <- list("NN" = FALSE,"RD_o" = FALSE,"RD" = FALSE,"SS_o" = FALSE,"SS" = TRUE,"pI" = FALSE,"OUR" = FALSE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE , "L2"=TRUE, "w_L2"= TRUE)

AUC_full_res<- list()


for(ind_data in 1:length(dataset)){
   
  # estimation of geo distance with methods TRUE in meth
  geo_estim = pairwise_geo_estimation(method=meth,
                                      true_data=NA,
                                      discrete_data=t(dataset[[ind_data]]$data),
                                      true_geo=NA,
                                      plot_true=FALSE,
                                      nb_proj=NA,
                                      grid=matrix(rep(dataset[[ind_data]]$reg_grid,length(dataset[[ind_data]]$reg_grid)),nrow=dataset[[ind_data]]$reg_grid,byrow = TRUE),
                                      reg_grid=dataset[[ind_data]]$reg_grid,
                                      common_grid_true=TRUE,
                                      Analytic_geo_available=FALSE,
                                      is_data_smoothed=TRUE)
  
  # Assess classificaton performance
  nr_train_test_splits = 500
  
  # Create balanced train-test partition taking into account the proportion of binary labels in the data
  # train.rows[,i] contains the training indices i-th split
  train.rows = createDataPartition(dataset[[ind_data]]$true_group, times = nr_train_test_splits, p = 0.53, list = FALSE)
  
  temp<- lapply(geo_estim,mean_auc,train.rows=train.rows,true_group=dataset[[ind_data]]$true_group,h=10)
  AUC_res <- matrix(unlist(temp,use.names=FALSE),ncol=length(geo_estim))
  boxplot(AUC_res,names=names(geo_estim),main=dataset[[ind_data]]$name,las=2,ylab="AUC")
  
  AUC_full_res[[dataset[[ind_data]]$name]]<-AUC_res
}

dev.off()


pdf("classification_results_grotwh_curve1.pdf",width=12,height=5)
par(mfrow=c(2,3))
for(i in 1:6){
  boxplot(AUC_full_res[[i]],las=2,ylab="AUC",names=c(substr(names(geo_estim)[1:4], 11, 17),substr(names(geo_estim)[5:6], 7, 10)),main=names(AUC_full_res)[i])
}
dev.off()

pdf("classification_results_grotwh_curve2.pdf",width=12,height=5)
par(mfrow=c(2,3))
for(i in 7:12){
  boxplot(AUC_full_res[[i]],las=2,ylab="AUC",names=c(substr(names(geo_estim)[1:4], 11, 17),substr(names(geo_estim)[5:6], 7, 10)),main=names(AUC_full_res)[i])
}
dev.off()




mean_auc <- function(train.rows,pairwise_distance_estimate,true_group,h){
  
  # train.rows created by createDataPartition
  # pairwise_distance_matrix samplesize by samplesize estimates of distance
  # true_group vector of length samplesize of binary class label
  nr_train_test_splits = dim(train.rows)[2]
  indices = 1:length(true_group)
  group1_indices = indices[true_group]
  aucs<-rep(0,nr_train_test_splits)
  class_err<-rep(0,nr_train_test_splits)
  
  for(split in 1:nr_train_test_splits){
    
    train_ind = train.rows[,split]
    
    # We classify the data using functional version of the Nadaraya–Watson kernel estimator for the class membership probabilities
    deno <- colSums(dnorm(h^(-1)*pairwise_distance_estimate[train_ind,-train_ind]))
    class1_count_hat<- colSums(dnorm(h^(-1)*pairwise_distance_estimate[train_ind[is.element(train_ind,group1_indices)],-train_ind]))
    class1_prob_hat<- class1_count_hat / deno
  
    predobj <- prediction(class1_prob_hat, true_group[-train_ind]*1)
    perf <- performance(predobj,"auc")
    aucs[split] = perf@y.values[[1]]
  }
  
  return(aucs)
}



