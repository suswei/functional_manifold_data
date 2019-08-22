choose_el<- function(list_in,ind1,ind2){
  list_in[[ind1]][ind2]
}

choose_el2<-function(ind_vec,list_in1){
  lapply(list_in1,choose_el,ind1=ind_vec[1],ind2=ind_vec[2])
}

create_boxplot<- function(sce_details,mat_ind,method_names,nr_methods){

  base_pattern = do.call(sprintf, c(list("sce=%d_SNR=%.01f_Kobs=%d"), sce_details))
  
  # TODO: fix hardcoding
  if(sce_details[2]==0.1){
    base_pattern_out = do.call(sprintf, c(list("sce=%d_SNR=low_Kobs=%d"), sce_details[1],sce_details[3]))
  }
  else{
    base_pattern_out = do.call(sprintf, c(list("sce=%d_SNR=high_Kobs=%d"), sce_details[1],sce_details[3]))
  }
  

  data_sce = list.files(pattern=paste("_", base_pattern,"_mc=*",sep=""))
  if (length(data_sce) == 0){
    return()
  }
  else{
    raw_results<-lapply(data_sce,readRDS)
    
    temp_res<-apply(mat_ind,1,choose_el2,list_in1=raw_results)
    temp_res2<- lapply(temp_res,unlist,use.names=FALSE)
    results<- matrix(unlist(temp_res2,use.names=FALSE),ncol=3*nr_methods)
    
    ### We plot one boxplot per method and per assesment measure
    pdf(paste(base_pattern_out,".pdf",sep=""),width=4,height=2)
    
    par(mfrow=c(1,3))
    # par(oma=c(3,3,3,3))
    boxplot(results[,1:nr_methods],main="MSE",names=method_names)
    boxplot(results[,(1+nr_methods):(2*nr_methods)],main="isometry",names=method_names)
    boxplot(results[,(2*nr_methods+1):(3*nr_methods)],main="Pearson",names=method_names)
    
    # mtext(base_pattern,outer=TRUE,cex=1.5)
    dev.off()
  }
}


setwd("./sim_results")

method_names=c("SS","s=2","s=3","L2")
nr_methods = length(method_names)

combination_res_assess<- cbind(rep(1:nr_methods,3),rep(1:3,each=nr_methods)) # nr_methods*3 combinations of methods and assesment measures
# TODO: fix hardcoding
combination_para <- expand.grid(c(5,2,4),c(0.5),c(30)) # combinations of the parameters, copied from compare_methods.R

apply(combination_para,1,create_boxplot,mat_ind=combination_res_assess,method_names=method_names,nr_methods=nr_methods)

setwd("../")

