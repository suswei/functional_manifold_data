#setwd('./Users/UQAM/Dropbox/Marie-Moi/Susan_project/manifold learning/Geo_dist_calculation/Cluster_results/spartanSim_clusterOutput')

choose_el<- function(list_in,ind1,ind2){
  list_in[[ind1]][ind2]
}

choose_el2<-function(ind_vec,list_in1){
  lapply(list_in1,choose_el,ind1=ind_vec[1],ind2=ind_vec[2])
}

create_boxplot<- function(sce_details,mat_ind,method_names,nr_methods){
  
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
  results<- matrix(unlist(temp_res2,use.names=FALSE),ncol=3*nr_methods)
  
  ### We plot one boxplot per method and per assesment measure
  pdf(paste("sce=",sce_details[1],
            "_samplesize=",sce_details[2],
            "_SNR=",sce_details[3],
            "_reg_sampling=",sce_details[4],
            ".pdf",sep=""),width=15,height=5)
  
  par(mfrow=c(1,3))
  par(oma=c(3,3,3,3))

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


# TODO: this needs to be passed into the script from compare_methods.R?
method_names=c("SS","Our3S1","Our3S2","Our3S3","L2")

nr_methods = length(method_names)


# TODO: a bit annoying that these are hardcoded, have to customise them according to compare_methods.R
combination_res_assess<- cbind(rep(1:nr_methods,3),rep(1:3,each=nr_methods)) # nr_methods*3 combinations of methods and assesment measures
combination_para <- expand.grid(c(5,2,4),c(0.1,0.5),c(100,30)) # combinations of the parameters

apply(combination_para,1,create_boxplot,mat_ind=combination_res_assess,method_names=method_names,nr_methods=nr_methods)



