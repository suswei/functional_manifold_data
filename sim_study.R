## Input 
# name 
# samplesize
# noise = sd of the noise
# reg_sampling : TRUE = regular sampling, FALSE = random uniform
# num_neigh = initial number of auto_isomap?

## Output
# plot of the true manlifold, of the observed one, and of the error in the geodesic estimation 
# for different choices of h



sim_study<- function(name,samplesize,noise,reg_sampling){
  
  par(mfrow=c(1,3))
  # load true and noisy manifold data from EuclideanExamples
  obj = EuclideanExamples(name, samplesize, noise, plotTrue=TRUE,reg_sampling,300)
  data = data.matrix(obj$data)
  true_mani  = data.matrix(obj$true_mani)
  true_geo = data.matrix(obj$true_geo)
  true_geo_tri_inf<- true_geo[lower.tri(true_geo, diag = FALSE)]
  
  
  # manifold estimation via subspace constrained mean shift (SCMS)
  scms_h=seq(0.1,3,by=0.1)
  Error_denoised_h= rep(0,length(scms_h))
  for(i in 1:length(scms_h)){
    denoised = scms$scms(data, scms_h[i])  #has same shape as data
    IsomapGdist_denoised = pyIsoAuto$get_auto_isomap_gdist(denoised)
    geo_dist<- IsomapGdist_denoised[lower.tri(IsomapGdist_denoised, diag = FALSE)]
    Error_denoised_h[i]=sqrt(sum((geo_dist -true_geo_tri_inf )^2))
  }
  
  # Calculus geo on the raw data
  IsomapGdist = pyIsoAuto$get_auto_isomap_gdist(data)
  geo_dist<-IsomapGdist[lower.tri(IsomapGdist, diag = FALSE)]
  err_raw=sqrt(sum((geo_dist -true_geo_tri_inf )^2))
  
  plot(scms_h,Error_denoised_h,ylim=c(min(Error_denoised_h,err_raw),max(Error_denoised_h,err_raw)),xlab='h',ylab='error',main=paste('n=',samplesize,' sampl.=',reg_sampling,' noise_sd=',noise,sep=''))
  abline(h=err_raw,col='red')
  
}







