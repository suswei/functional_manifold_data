#setwd('./Users/UQAM/Dropbox/Marie-Moi/Susan_project/manifold learning/Geo_dist_calculation/Cluster_results/spartanSim_clusterOutput')

choose_el<- function(list_in,ind1,ind2){
  list_in[[ind1]][ind2]
}

choose_el2<-function(ind_vec,list_in1){
  lapply(list_in1,choose_el,ind1=ind_vec[1],ind2=ind_vec[2])
}

create_boxplot<- function(sce_details,mat_ind){
  
  data_sce = list.files(pattern=paste("_sce=",sce_details[1],
                                      "_samplesize=",sce_details[2],
                                      "_SNR=",sce_details[3],
                                      "_reg_sampling=",sce_details[4]
                                      ,"_mc=*",sep=""))
  if (length(data_sce) == 0){
    return()
  }
  else{
  raw_results<-lapply(data_sce,readRDS)
  
  temp_res<-apply(mat_ind,1,choose_el2,list_in1=raw_results)
  temp_res2<- lapply(temp_res,unlist,use.names=FALSE)
  results<- matrix(unlist(temp_res2,use.names=FALSE),ncol=18)
  
  ### We plot one boxplot per method and per assesment measure
  pdf(paste("sce=",sce_details[1],
            "_samplesize=",sce_details[2],
            "_SNR=",sce_details[3],
            "_reg_sampling=",sce_details[4],
            ".pdf",sep=""),width=15,height=5)
  
  par(mfrow=c(1,3))
  par(oma=c(2,2,2,2))

  meth <- list("NN" = TRUE,"RD_o" = TRUE,"RD" = TRUE,"SS_o" = TRUE,"SS" = TRUE,"pI" = FALSE,"OUR" = TRUE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE )# see pairwise_geo_estimation for more info
  
  # TODO: this needs to be passed into the script from compare_methods.R?
  method_names=c("NN","RD_o", "RD", "SS_o","SS","OurS1","OurS2","OurS3","Our3_S1","Our3_S2","Our3_S3")
  nr_methods = length(method_names)
  
  # TODO: VERY annoying that these are hardcoded, have to customise them according to compare_methods.
  boxplot(results[,1:nr_methods],main="relative MSE",names=method_names)
  boxplot(results[,(1+nr_methods):(2*nr_methods)],main="AUC of entrywise epsilon-isometry",names=method_names)
  boxplot(results[,(2*nr_methods+1):(3*nr_methods)],main="Pearson correlation",names=method_names)
  mtext(paste("sce=",sce_details[1],"_samplesize=",sce_details[2],
              "_SNR=",sce_details[3],
              "_reg_sampling=",sce_details[4],sep=""),outer=TRUE,cex=1.5)
  dev.off()
  }
}

# TODO: a bit annoying that these are hardcoded, have to customise them according to compare_methods.R
combination_res_assess<- cbind(rep(1:nr_methods,3),rep(1:3,each=nr_methods)) # 7*3 combinations of methods and assesment measures
combination_para <- expand.grid(c(1,2,4),c(100,250),c(0.1,0.5),c(0,1)) # combinations of the parameters

apply(combination_para,1,create_boxplot,mat_ind=combination_res_assess)



