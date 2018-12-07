# TODO Marie, this function should call assess_goodness_estimation!
##### Estimation of geodesic distances

# We estimate the geodesic matrix obtained via
# a) Floyd's Algorithm on true data (just to make sure that the theoretical geo. is calculated correctly)
# b) Floyd's Algorithm on raw data
# c) Floyd's Algorithm on smooth data (only for FD)
# d) p-isomap of Muller on smooth data
# e) mds of smooth data + scsm + Floyd's Algorithm (only for FD) 
# f) scsm + Floyd's Algorithm (only for ED)
# g) random proj. of smooth data + scsm + Floyd's Algorithm (only for FD)


## Note that a) and b) can only be used if the data are observed on a common grid
## Note that c) can only be used for functional data
## In the case of Eucledian data, d), e) and f) use directly the raw data

## Input
# true_data : samplesize x K matrix containing the original data (no noise)
# discrete_data : samplesize x K matrix containing the observed data
# true_geo : samplesize x samplesize matrix containing true pairwise geo disctance
# s : dimension use for mds and random projections
# plot_true : if TRUE plot of geo estimation for each method
# FD_true : if TRUE functional data, otherwise euclidean data
# grid : samplesize x K matrix containing the grid on which each data is observed (required only if FD_true =1)
# common_grid_true : if 1 grid is common for every curve, if 0 different grid for every curve

## Output
# list of estimated geodesic distance :
# estim_geo_true_data = method a)
# estim_geo_noisy_data = method b)
# estim_geo_smooth_data = method c)
# estim_geo_penalized_isomap = method d)
# estim_geo_mds_scms = method e)
# estim_geo_scms = method f)
# estim_geo_RP_scms = method g)


#Example of call for FD
# Estim<- Geo_estimation(data$true_data,data$discrete_data,data$true_geo,TRUE,TRUE,3,20,data$grid,data$reg_grid,1)
#Example of call for FD
# Estim<- Geo_estimation(data$true_data,data$discrete_data,data$true_geo,TRUE,FALSE)


Geo_estimation <- function(true_data,discrete_data,true_geo,plot_true,FD_true,s,nb_proj,grid,reg_grid,common_grid_true=0){
  
  samplesize = nrow(true_data)
  K = ncol(true_data)
  norm_true_geo = sqrt(sum(true_geo^2))
  list1_to_return=list()
  list2_to_return=list()
  list3_to_return=list()
  list4_to_return=list()
  
  ## Apply method a) and b) for FD observed on a common grid or Euclidean data
  if(!FD_true||common_grid_true==1){ 
    
    ### a) Estimation from the true data 
    
    print("method a running")
    
    if(FD_true){
      a = reg_grid[1]
      b = reg_grid[K]
      true_data_tmp = (sqrt((b-a)/K))*true_data
    } else{
      true_data_tmp = true_data
    }
    
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(true_data_tmp)
    num_neigh_true=seq(num_neigh_min,samplesize/2,by=2)
    
    # Calculate the geo matrix for each number of neighbors and keep the one that gives the minimal error
    Error_true_mani_K= rep(0,length(num_neigh_true))
    for(j in 1:length(num_neigh_true)){
      IsomapGdist = pyIso$getIsomapGdist(true_data_tmp,num_neigh_true[j])
      Error_true_mani_K[j]=sqrt(sum((IsomapGdist -true_geo )^2))/norm_true_geo
    }
    ind_op_true=min(which(Error_true_mani_K==min(Error_true_mani_K)))
    estim_geo_true_data = pyIso$getIsomapGdist(true_data_tmp,num_neigh_true[ind_op_true])
    
    ### b) Direct estimation from the noisy data
    
    print("method b running")
    
    if(FD_true){
      a = reg_grid[1]
      b = reg_grid[K]
      discrete_data_tmp = (sqrt((b-a)/K))*discrete_data
    } else {
      discrete_data_tmp = discrete_data
    }
    
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(discrete_data_tmp)
    num_neigh_noisy=seq(num_neigh_min,samplesize/2,by=2)
    
    # Calculate the geo matrix for each number of neighbors and keep the one that gives the minimal error
    Error_noisy_mani_K= rep(0,length(num_neigh_noisy))
    for(j in 1:length(num_neigh_noisy)){
      IsomapGdist = pyIso$getIsomapGdist(discrete_data_tmp,num_neigh_noisy[j])
      Error_noisy_mani_K[j]=sqrt(sum((IsomapGdist -true_geo )^2))/norm_true_geo
    }
    ind_op_noisy=min(which(Error_noisy_mani_K==min(Error_noisy_mani_K)))
    estim_geo_noisy_data = pyIso$getIsomapGdist(discrete_data_tmp,num_neigh_noisy[ind_op_noisy])
    
    if(plot_true){
      par(mfrow=c(1,2))
      image.plot(estim_geo_true_data,main=paste('Err true data = ',signif(min(Error_true_mani_K),digits =4),', nb nei. =',num_neigh_true[ind_op_true],sep=''))
      image.plot(estim_geo_noisy_data,main=paste('Err noisy data = ',signif(min(Error_noisy_mani_K),digits =4),', nb nei. =',num_neigh_noisy[ind_op_noisy],sep=''))
    }
    
    # Update list of estimation geo
    list1_to_return = list("estim_geo_true_data"=estim_geo_true_data,"estim_geo_noisy_data"=estim_geo_noisy_data)
    
  }
  
  
  ## Apply method c) for functional data 
  if(FD_true){
    
    print("method c running")
    
    # Smooth the observed data using b-spline and evaluate the full curve on the same regular grid
    smooth_data = smooth_FD_bspline(discrete_data,grid,reg_grid)
    a = reg_grid[1]
    b = reg_grid[K]
    smooth_data_tmp = (sqrt((b-a)/K))*smooth_data
    ### c) Estimation from the smooth curve
    
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(smooth_data_tmp)
    num_neigh_smooth=seq(num_neigh_min,samplesize/2,by=2)
    
    # Calculate the geo matrix for each number of neighbors and keep the one that gives the minimal error
    Error_smooth_mani_K= rep(0,length(num_neigh_smooth))
    for(j in 1:length(num_neigh_smooth)){
      IsomapGdist = pyIso$getIsomapGdist(smooth_data_tmp,num_neigh_smooth[j])
      Error_smooth_mani_K[j]=sqrt(sum((IsomapGdist -true_geo )^2))/norm_true_geo
    }
    ind_op_smooth=min(which(Error_smooth_mani_K==min(Error_smooth_mani_K)))
    estim_geo_smooth_data = pyIso$getIsomapGdist(smooth_data_tmp,num_neigh_smooth[ind_op_smooth])
    
    if(plot_true){
      par(mfrow=c(1,2))
      matplot(reg_grid,t(smooth_data),main="Smoothed data",type='l', col=rainbow(samplesize))
      image.plot(estim_geo_smooth_data,main=paste('Err smoothed data = ',signif(min(Error_smooth_mani_K),digits = 4),'nb nei. =',num_neigh_smooth[ind_op_smooth],sep=''))
    }
    
    ## change name smooth_data to discrete_data (and "normalize") to use the same procedure 
    ## as for Euclidean data for method d), e) and f)
    discrete_data = smooth_data_tmp
    
    # Update list of estimation geo
    list2_to_return = list("estim_geo_smooth_data"=estim_geo_smooth_data)
  }
  
  ### d) p-isomap of Muller
  
  print("method d running")
  
  # Find a grid of possible values for the number of neigbors
  num_neigh_min=py_min_neigh$get_min_num_neighbors(discrete_data)
  num_neigh=seq(num_neigh_min,samplesize/2,by=2)
  
  #Define candidate values for delta
  delta_can=c(0,0.02,0.05,0.1,0.2)
  
  write.table(delta_can,file='./possible_delta.txt',col.names = FALSE, row.names = FALSE)
  write.table(discrete_data,file='./discrete_data.txt',col.names = FALSE, row.names = FALSE)
  write.table(num_neigh,file="./possible_K.txt",col.names = FALSE, row.names = FALSE)
  write.table(c(samplesize,K),file="./n_K.txt",col.names = FALSE, row.names = FALSE)
  write.table(true_geo,file="./true_geo.txt",col.names = FALSE, row.names = FALSE)
  
  # Run penalized-isomap.m in matlab
  run_matlab_script("./run_isomap_p.m")
  
  temp=as.vector(read.table(file='./manidis.txt',header=FALSE))
  estim_geo_penalized_isomap=matrix(temp[,1],ncol=samplesize)
  temp2=read.table(file='./err_delta_neigh.txt',header=FALSE)
  err_delta_neigh=matrix(temp2[,1],ncol=length(num_neigh))
  delta_min_tmp=which(err_delta_neigh==min(err_delta_neigh),arr.ind = TRUE)
  delta_min=delta_can[delta_min_tmp[1,1]]
  
  ### e) Functional data : use mds to obtain a s-d version of each data point, use SCMS to estimate the manifold 
  ###    and estimation geo with Floyd's algo.
  ### f) Euclidean data : use SCMS to estimate the manifold 
  ###    and estimation geo with Floyd's algo.
  print("method e/f running")
  
  if(FD_true){
    # Projection of the data in dim s using mds
    euc_dis_mat=dist(discrete_data)
    mds_mat = cmdscale(euc_dis_mat,eig=TRUE, k=s) 
    projected = mds_mat$points
  } else{
    projected = discrete_data
  }
  
  scms_h=seq(0.1,2,by=0.2)
  Error_mds_h= matrix(0,nrow=length(scms_h),ncol=2)
  for(bw in 1:length(scms_h)){
    # use scms to estimate a manifold in dim s
    denoised = scms$scms(projected, scms_h[bw])
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(denoised)
    num_neigh=seq(num_neigh_min,samplesize/2,by=2)
    Error_mds_K= rep(0,length(num_neigh))
    for(j in 1:length(num_neigh)){
      IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,num_neigh[j])
      Error_mds_K[j]=sqrt(sum((IsomapGdist_denoised -true_geo)^2))/norm_true_geo
    }
    ind_neigh=min(which(Error_mds_K==min(Error_mds_K)))
    Error_mds_h[bw,]=c(Error_mds_K[ind_neigh],num_neigh[ind_neigh])
  }
  ind=min(which(Error_mds_h[,1]==min(Error_mds_h[,1])))
  denoised = scms$scms(projected, scms_h[ind])
  estim_geo_mds_scms = pyIso$getIsomapGdist(denoised,Error_mds_h[ind,2])
  
  if(FD_true){
    list3_to_return = list("estim_geo_penalized_isomap"=estim_geo_penalized_isomap,"pisomap_delta" = delta_min,"estim_geo_mds_scms"=estim_geo_mds_scms)
  } else{
    list3_to_return = list("estim_geo_penalized_isomap"=estim_geo_penalized_isomap,"estim_geo_scms"=estim_geo_mds_scms)
  }
  
  if(plot_true){
    par(mfrow=c(1,2))
    image.plot(estim_geo_penalized_isomap,main=paste('err p-isomap = ',signif(min(err_delta_neigh),digits =4),', delta = ',delta_min,sep=''))
    if(FD_true){
      image.plot(estim_geo_mds_scms,main=paste('err mds + scms = ',signif(min(Error_mds_h[,1]),digits = 4),sep=''))
    }else{
      image.plot(estim_geo_mds_scms,main=paste('err scms = ',signif(min(Error_mds_h[,1]),digits = 4),sep=''))
    }
  }
  
  ### f) use random projection to obtain a s-d version of each data point, use SCMS to estimate the manifold 
  ###    and estimation geo with Floyd's algo.
  if(FD_true){
    print("method f running")
  
    scms_h=seq(0.1,2,by=0.2)
    geo_proj<-matrix(0,nrow=nb_proj,ncol=samplesize*(samplesize-1)/2)
    op_h_proj<-rep(0,nb_proj)
    op_num_nei_proj<-rep(0,nb_proj)
    err_proj<-rep(0,nb_proj)
  
    for(proj in 1:nb_proj){
      # use random projection to represent data in dim s
      randprojmat = matrix(rnorm(s*K,mean=0,sd=1/sqrt(s)), nrow=s, ncol=K) 
      projected = randprojmat %*% t(discrete_data)
      Error_RP_h= matrix(0,nrow=length(scms_h),ncol=2)
      for(bw in 1:length(scms_h)){
        # use scms to estimate a manifold in dim s
        denoised = scms$scms(t(projected), scms_h[bw])
        # Find a grid of possible values for the number of neigbors
        num_neigh_min=py_min_neigh$get_min_num_neighbors(denoised)
        num_neigh=seq(num_neigh_min,samplesize/2,by=2)
        Error_RP_K= rep(0,length(num_neigh))
        for(j in 1:length(num_neigh)){
          IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,num_neigh[j])
          Error_RP_K[j]=sqrt(sum((IsomapGdist_denoised -true_geo)^2))/norm_true_geo
        }
        ind_neigh=min(which(Error_RP_K==min(Error_RP_K)))
        Error_RP_h[bw,]=c(Error_RP_K[ind_neigh],num_neigh[ind_neigh])
      }
      ind=min(which(Error_RP_h[,1]==min(Error_RP_h[,1])))
      denoised = scms$scms(t(projected), scms_h[ind])
      IsomapGdist_denoised = pyIso$getIsomapGdist(denoised,Error_RP_h[ind,2])
      geo_proj[proj,]<-IsomapGdist_denoised[lower.tri(IsomapGdist_denoised, diag = FALSE)]
      op_h_proj[proj]=scms_h[ind]
      op_num_nei_proj[proj]=Error_RP_h[ind,2]
      err_proj[proj]=sqrt(sum((IsomapGdist_denoised -true_geo )^2))/norm_true_geo
    }
  
    ensemble_geo<-apply(geo_proj,2,median)
    ens_geo_mat<-matrix(0,samplesize,samplesize)
    ens_geo_mat[lower.tri(ens_geo_mat, diag = FALSE)]=ensemble_geo
    estim_geo_RP_scms= ens_geo_mat + t(ens_geo_mat)
    err_RP<-sqrt(sum((estim_geo_RP_scms -true_geo )^2))/norm_true_geo
  
    if(plot_true){
      par(mfrow=c(1,1))
      image.plot(estim_geo_RP_scms,main=paste('err RP + scms = ',signif(err_RP,digits = 4),sep=''))
    }
    list4_to_return = list("estim_geo_RP_scms"=estim_geo_RP_scms)
  }
  
  return(c(list1_to_return,list2_to_return,list3_to_return,list4_to_return))
}