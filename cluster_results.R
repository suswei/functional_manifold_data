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
    
    pdf(paste(base_pattern_out,"_MSE",".pdf",sep=""),width=16/4,height=3)
    par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.5,1.5,1.5,1.5))
    # par(
    #   # Change the colors
    #   col.main="red", col.lab="blue", col.sub="black",
    #   # Titles in italic and bold
    #   font.main=4, font.lab=4, font.sub=4,
    #   # Change font size
    #   cex.main=0.5, cex.lab=0.5, cex.sub=0.5
    # )
    boxplot(results[,1:nr_methods],main="MSE",names=method_names,cex.main=0.75,cex.names=0.75,cex.axis=0.75,cex.lab=0.75)
    dev.off()
    
    pdf(paste(base_pattern_out,"_isometry",".pdf",sep=""),width=16/4,height=3)
    par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.5,1.5,1.5,1.5))
    boxplot(results[,(1+nr_methods):(2*nr_methods)],main="Isometry",names=method_names,cex.main=0.75,cex.names=0.75,cex.axis=0.75,cex.lab=0.75)
    dev.off()
    
    pdf(paste(base_pattern_out,"_Pearson",".pdf",sep="") ,width=16/4,height=3)
    par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.5,1.5,1.5,1.5))
    boxplot(results[,(2*nr_methods+1):(3*nr_methods)],main="Pearson",names=method_names,cex.main=0.75,cex.names=0.75,cex.axis=0.75,cex.lab=0.75)
    dev.off()
    
    # mtext(base_pattern,outer=TRUE,cex=1.5)
  }
}
