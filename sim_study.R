# Function to generate data
## Input 
# name 
# samplesize
# noise = sd of the noise
# reg_sampling : TRUE = regular sampling, FALSE = random uniform
## Output
# 3 files .txt : one with the true manifold, one with the observed one (with noise), and the true geodesic distance matrix
##### YOU NEED TO CREATE A FOLDER IN THE CURRECNT DIRRECTORY CALLED simulated_data BEFORE USING THE FUNCTION
generate_data<-function(name,samplesize,SNR,reg_sampling){
  obj = EuclideanExamples(name, samplesize, SNR, plotTrue=FALSE,reg_sampling,300)
  data = data.matrix(obj$data)
  true_mani  = data.matrix(obj$true_mani)
  true_geo = data.matrix(obj$true_geo)
  write.table(data,file=paste("./simulated_data/",name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'_data.txt',sep=''),col.names = FALSE,row.names = FALSE)
  write.table(true_mani,file=paste("./simulated_data/",name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'_mani.txt',sep=''),col.names = FALSE,row.names = FALSE)
  write.table(true_geo,file=paste("./simulated_data/",name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'_geo.txt',sep=''),col.names = FALSE,row.names = FALSE)
}


#Function to analyse the data
## Input 
# name 
# samplesize
# noise = sd of the noise
# reg_sampling : TRUE = regular sampling, FALSE = random uniform
# do_plot : TRUE or FALSE
# do_pdf : TRUE or FALSE
## Output
# if plot = TRUE : plot of the true manlifold, of the observed one, and of the error in the geodesic estimation for different choices of h
# if pdf=TRUE same plot in a pdf file

sim_study<- function(name,samplesize,SNR,reg_sampling,do_plot,do_pdf){
  
  # load true and noisy manifold data from the folder Simulated_data
  data<- as.matrix(read.table(file=paste("./simulated_data/",name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'_data.txt',sep='')))
  true_mani<- as.matrix(read.table(file=paste("./simulated_data/",name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'_mani.txt',sep='')))
  true_geo<- as.matrix(read.table(file=paste("./simulated_data/",name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'_geo.txt',sep='')))
  true_geo_tri_inf<- true_geo[lower.tri(true_geo, diag = FALSE)]
  
  
  #manifold estimation via subspace constrained mean shift (SCMS)
  scms_h=seq(0.1,1.5,by=0.1)
  Error_denoised_h= rep(0,length(scms_h))
  num_neigh=seq(2,15,by=1)
  
  for(i in 1:length(scms_h)){
    Error_denoised_K= rep(0,length(num_neigh))
    denoised = scms$scms(data, scms_h[i])  #has same shape as data
    for(j in 1:length(num_neigh)){
      IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,num_neigh[j])
      geo_dist<- IsomapGdist_denoised[lower.tri(IsomapGdist_denoised, diag = FALSE)]
      Error_denoised_K[j]=sqrt(sum((geo_dist -true_geo_tri_inf )^2))
    }
    Error_denoised_h[i]=min(Error_denoised_K)
  }
  
  
  # Calculus geo on the raw data
  Error_denoised_K= rep(0,length(num_neigh))
  for(j in 1:length(num_neigh)){
    IsomapGdist = pyIso$getIsomapGdist(data,num_neigh[j])
    geo_dist<- IsomapGdist[lower.tri(IsomapGdist, diag = FALSE)]
    Error_denoised_K[j]=sqrt(sum((geo_dist -true_geo_tri_inf )^2))
  }
  err_raw=min(Error_denoised_K)
  
  if(do_plot==TRUE){
    plot(scms_h,Error_denoised_h,type="b",ylim=c(min(Error_denoised_h,err_raw),max(Error_denoised_h,err_raw)),xlab='h',ylab='error',main=paste('n=',samplesize,' sampl.=',reg_sampling,' SNR=',SNR,sep=''))
    abline(h=err_raw,col='red')
  }
  if(do_pdf==TRUE){
    pdf(file=paste(name,'_n=',samplesize,'_SNR=',SNR,'_reg_sam=',reg_sampling,'.pdf',sep=''))
    plot(scms_h,Error_denoised_h,type="b",ylim=c(min(Error_denoised_h,err_raw),max(Error_denoised_h,err_raw)),xlab='h',ylab='error',main=paste('n=',samplesize,' sampl.=',reg_sampling,' SNR=',SNR,sep=''))
    abline(h=err_raw,col='red')
    dev.off()
  }
  
}
