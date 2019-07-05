##### Estimation of pairwise geodesic distances for functional data

# We estimate the geodesic matrix obtained via
# NN: Floyd's Algorithm on noiseless data (just to make sure that the theoretical geo. is calculated correctly)
# RD_o: Floyd's Algorithm on noisy data (oracle for num neigh)
# RD: Floyd's Algorithm on noisy data (no oracle)
# SS_o: Floyd's Algorithm on smooth data (oracle for num neigh)
# SS: Floyd's Algorithm on smooth data (no oracle)
# pI: p-isomap of Muller on smooth data
# OUR: mds of smooth data + scms + Floyd's Algorithm (oracle for h and num neigh)
# OUR2: mds of smooth data + scms + robust isomap (oracle for h)
# OUR3: mds of smooth data + scms h heuri + robust isomap (no oracle)
# RP : random proj. of smooth data + scms + Floyd's Algorithm 
# L2 : pairwise L2 distances on smooth data
# w_L2 : weigthed L2 defined in Chen et al (Biometrics)


# TODO: what if data are not observed on a common grid and NN, RD, RD_o are called? Are there warnings?
## Note that NN, RD, RD_o can only be used if the data are observed on a common grid

## Input
# method : list of dimension 8 indicating which estimation methods should be run. Warning : NN and RD only works if common_grid_true=1
# true_data [optional, bypass with NA]: samplesize x K matrix containing the original data (no noise)
# discrete_data : samplesize x K matrix containing the observed data
# true_geo [optional, bypass with NA]: samplesize x samplesize matrix containing true pairwise geo disctance
# s : dimension use for mds and random projections
# plot_true : if TRUE plot of geo estimation for each method
# nb_proj : number of projection used for RP method
# grid : samplesize x K matrix containing the grid on which each data is observed 
# reg_grid : vector of dim K containing a regular grid of K points on [a,b] 
# K_dense : size of the regular grid on which smoothed curve are evaluated 
# common_grid_true : if 1 grid is common for every curve, if 0 different grid for every curve
# Analytic_geo_available : indicate if we have access to the true pairwise geodesic distances
# is_data_smoothed : TRUE data has already been smoothed, otherwise false



## Output
# list of estimated geodesic distance :
# estim_geo_NN = method NN
# estim_geo_RD_o = method RD_o
# estim_geo_RD = method RD
# estim_geo_SS_o = method SS_o
# estim_geo_SS = method SS
# estim_geo_pI = method pI
# estim_geo_OUR = method OUR
# estim_geo_OUR2 = method OUR2
# estim_geo_OUR3 = method OUR3
# estim_geo_RP = method RP
# estim_L2 =  method L2 
# estim_w_L2 = method w_L2


# Example of call for functional data
# data<- sim_functional_data(sce=2,samplesize=100)
# meth <- list("NN" = TRUE,"RD_o" = TRUE,"RD" = TRUE,"SS_o" = TRUE,"SS" = TRUE,"pI" = FALSE,"OUR" = TRUE,"OUR2" = FALSE,"OUR3"=TRUE,"RP" = FALSE, "L2" = TRUE, "w_L2" = TRUE )
# Estim<- pairwise_geo_estimation(meth,data$noiseless_data,data$noisy_data,data$analytic_geo,TRUE,TRUE,20,data$grid,data$reg_grid,100,1,FALSE,FALSE)


pyIso = import_from_path("getIsomapGdist",path='.') # TODO: consider replacing python Isomap with R implementation of Floyd's algorithm? https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
py_min_neigh = import_from_path("get_min_num_neighbors",path='.')
scms = import_from_path("scms",path='.') # TODO: consider replacing python scms with R version, see https://sites.google.com/site/yenchicr/algorithm

pairwise_geo_estimation <- function(method,
                                    true_data,
                                    discrete_data,
                                    true_geo,
                                    plot_true,
                                    nb_proj=NULL,
                                    grid,
                                    reg_grid,
                                    K_dense,
                                    common_grid_true,
                                    Analytic_geo_available=TRUE,
                                    is_data_smoothed = FALSE){
  
  samplesize = nrow(discrete_data)
  K = ncol(discrete_data)
  a = reg_grid[1]
  b = reg_grid[K]
  reg_grid_dense<- seq(a,b,length.out=K_dense)
  list_to_return<- list()

  
  if(is_data_smoothed) {
    smooth_data = discrete_data
    smooth_data_K <- discrete_data
    # Normalise the data
    smooth_data <- (sqrt((b-a)/(K-1)))*smooth_data
    smooth_data_K <- (sqrt((b-a)/(K-1)))*smooth_data
  } else {
    smooth_fd<- smooth_FD_bspline(discrete_data,grid,reg_grid)
    smooth_data <- t(eval.fd(reg_grid_dense,smooth_fd))
    smooth_data_K<- t(eval.fd(reg_grid,smooth_fd))
    # Normalise the data
    smooth_data <- (sqrt((b-a)/(K_dense-1)))*smooth_data
    smooth_data_K <- (sqrt((b-a)/(K-1)))*smooth_data_K
  }
  
  # Normalise the data
  discrete_data <- (sqrt((b-a)/(K-1)))*discrete_data
  
  # plot noisy and smooth data
  if(plot_true){
      matplot(reg_grid,t(discrete_data),main="Raw data",type='l', col=rainbow(samplesize))
      matplot(reg_grid,t(smooth_data_K),main="Smoothed data",type='l', col=rainbow(samplesize))
  }
  
  
  if(Analytic_geo_available){
    norm_true_geo = sqrt(sum(true_geo^2))
  }
  
  if(method$NN){
    
    ### Estimation from noiseless data
    
    print("Floyd's Algorithm on noiseless data")
    true_data <- (sqrt((b-a)/(K-1)))*true_data
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(true_data)
    num_neigh_true=seq(num_neigh_min,samplesize/2,by=2)
    
    # Calculate the geo matrix for each number of neighbors and keep the one that gives the minimal error
    Error_true_mani_K= rep(0,length(num_neigh_true))
    for(j in 1:length(num_neigh_true)){
      IsomapGdist = pyIso$getIsomapGdist(true_data,num_neigh_true[j])
      Error_true_mani_K[j]=sqrt(sum((IsomapGdist -true_geo )^2))/norm_true_geo
    }
    ind_op_true=min(which(Error_true_mani_K==min(Error_true_mani_K)))
    estim_geo_true_data = pyIso$getIsomapGdist(true_data,num_neigh_true[ind_op_true])
    
    list_to_return[["estim_geo_NN"]] <- estim_geo_true_data
    
    if(plot_true){
      image.plot(estim_geo_true_data,main=paste('Err noiseless data = ',signif(min(Error_true_mani_K),digits =4),', nb nei. =',num_neigh_true[ind_op_true],sep=''))
    }
  }
  
  if(method$RD_o){
    ### Direct estimation from the noisy data
    
    print("Floyd's Algorithm on raw data")
    
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(discrete_data)
    num_neigh_noisy=seq(num_neigh_min,samplesize/2,by=2)
    
    # Calculate the geo matrix for each number of neighbors and keep the one that gives the minimal error
    Error_noisy_mani_K= rep(0,length(num_neigh_noisy))
    for(j in 1:length(num_neigh_noisy)){
      IsomapGdist = pyIso$getIsomapGdist(discrete_data,num_neigh_noisy[j])
      Error_noisy_mani_K[j]=sqrt(sum((IsomapGdist -true_geo )^2))/norm_true_geo
    }
    ind_op_noisy=min(which(Error_noisy_mani_K==min(Error_noisy_mani_K)))
    estim_geo_noisy_data = pyIso$getIsomapGdist(discrete_data,num_neigh_noisy[ind_op_noisy])
    
    if(plot_true){
      image.plot(estim_geo_noisy_data,main=paste('Err noisy data = ',signif(min(Error_noisy_mani_K),digits =4),', nb nei. =',num_neigh_noisy[ind_op_noisy],sep=''))
    }
    
    # Update list of estimation geo
    list_to_return[["estim_geo_RD_o"]] <- estim_geo_noisy_data
    
  }
  
  if(method$RD){
    ### Direct estimation from the noisy data
    
    print("Floyd's Algorithm on raw data (no oracle)")
    estim_geo_noisy_data<-robust_iso(discrete_data)
    
    if(plot_true){
      err<- sqrt(sum((estim_geo_noisy_data -true_geo)^2))/norm_true_geo
      image.plot(estim_geo_noisy_data,main=paste('Err noisy data = ',signif(err,digits =4),sep=''))
    }
    
    # Update list of estimation geo
    list_to_return[["estim_geo_RD"]] <- estim_geo_noisy_data
    
  }
  

  if(method$SS_o){
    ## Estimation from smooth data
    
    print("Floyd's Algorithm on smooth data")
   
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(smooth_data)
    num_neigh_smooth=seq(num_neigh_min,samplesize/2,by=2)
    
    # Calculate the geo matrix for each number of neighbors and keep the one that gives the minimal error
    Error_smooth_mani_K= rep(0,length(num_neigh_smooth))
    for(j in 1:length(num_neigh_smooth)){
      IsomapGdist = pyIso$getIsomapGdist(smooth_data,num_neigh_smooth[j])
      Error_smooth_mani_K[j]=sqrt(sum((IsomapGdist -true_geo )^2))/norm_true_geo
    }
    ind_op_smooth=min(which(Error_smooth_mani_K==min(Error_smooth_mani_K)))
    estim_geo_smooth_data = pyIso$getIsomapGdist(smooth_data,num_neigh_smooth[ind_op_smooth])
    
    if(plot_true){
      image.plot(estim_geo_smooth_data,main=paste('Err smoothed data = ',signif(min(Error_smooth_mani_K),digits = 4),'nb nei. =',num_neigh_smooth[ind_op_smooth],sep=''))
    }
    
    # Update list of estimation geo
    list_to_return[["estim_geo_SS_o"]] <- estim_geo_smooth_data
  }
  
  if(method$SS){
    ## Estimation from smooth data
    
    print("Floyd's Algorithm on smooth data (no oracle)")
    estim_geo_smooth_data = robust_iso(smooth_data)
    
    if(plot_true){
      err<- sqrt(sum((estim_geo_smooth_data -true_geo )^2))/norm_true_geo
      image.plot(estim_geo_smooth_data,main=paste('Err smoothed data = ',signif(err,digits = 4),sep=''))
    }
    
    # Update list of estimation geo
    list_to_return[["estim_geo_SS"]] <- estim_geo_smooth_data
  }
  
  if(method$pI){
    ### p-isomap of Muller
  
    print("p-isomap of Muller on smooth data")
  
    # Find a grid of possible values for the number of neigbors
    num_neigh_min=py_min_neigh$get_min_num_neighbors(smooth_data)
    num_neigh=seq(num_neigh_min,samplesize/2,by=2)
  
    #Define candidate values for delta
    delta_can=c(0,0.02,0.05,0.1,0.2)
  
    write.table(delta_can,file='./possible_delta.txt',col.names = FALSE, row.names = FALSE)
    write.table(smooth_data,file='./discrete_data.txt',col.names = FALSE, row.names = FALSE)
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
  
    list_to_return[["estim_geo_pI"]]<- estim_geo_penalized_isomap
    
    if(plot_true){
      image.plot(estim_geo_penalized_isomap,main=paste('err p-isomap = ',signif(min(err_delta_neigh),digits =4),', delta = ',delta_min,sep=''))
    }
  }

  
  if(method$OUR){
  ###   use mds to obtain a s-dim version of each data point, use SCMS to estimate the manifold
  ###   and estimation geo with Floyd's algo.

    print("mds of smooth data + scms + Floyd's Algorithm (oracle for num. neigh and h)")

    err_s<- rep(0,3)
    scms_h=seq(0.05,1.5,by=0.1)
    for(s in 1:3){
      # Projection of the data in dim s using mds
      euc_dis_mat=dist(smooth_data)
      mds_mat = cmdscale(euc_dis_mat,eig=TRUE, k=s)
      projected = mds_mat$points

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
      list_to_return[[paste("estim_geo_OUR_s",s,sep='')]] = pyIso$getIsomapGdist(denoised,Error_mds_h[ind,2])
      err_s[s]<-min(Error_mds_h[,1])
    }

    if(plot_true){
      par(mfrow=c(1,3))
      image.plot(list_to_return$estim_geo_OUR_s1,main=paste('OUR  (s=1) = ',signif(err_s[1],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_OUR_s2,main=paste('OUR (s=2) = ',signif(err_s[2],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_OUR_s3,main=paste('OUR (s=3) = ',signif(err_s[3],digits = 4),sep=''))
    }
  }
  
  
  if(method$OUR2){
    ###   use mds to obtain a s-dim version of each data point, use SCMS to estimate the manifold
    ###   with oracle h and estimation geo with Floyd's algo.
    
    print("mds of smooth data + scms + Robust Iso (oracle h)")
    
      err_s<- rep(0,3)
      scms_h=seq(0.05,1.5,by=0.1)
      for(s in 1:3){
        # Projection of the data in dim s using mds
        euc_dis_mat=dist(smooth_data)
        mds_mat = cmdscale(euc_dis_mat,eig=TRUE, k=s)
        projected = mds_mat$points
        
        Error_mds_h= rep(0,nrow=length(scms_h))
        for(bw in 1:length(scms_h)){
          # use scms to estimate a manifold in dim s
          denoised = scms$scms(projected, scms_h[bw])
          estim_bw<-robust_iso(denoised)
          Error_mds_h[bw]=sqrt(sum((estim_bw -true_geo)^2))/norm_true_geo
        }
        ind=min(which.min(Error_mds_h))
        denoised = scms$scms(projected, scms_h[ind])
        list_to_return[[paste("estim_geo_OUR2_s",s,sep='')]] = robust_iso(denoised)
        err_s[s]<-min(Error_mds_h)
      }
    
    if(plot_true){
      par(mfrow=c(1,3))
      image.plot(list_to_return$estim_geo_OUR2_s1,main=paste('OUR2 (s=1) = ',signif(err_s[1],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_OUR2_s2,main=paste('OUR2 (s=2) = ',signif(err_s[2],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_OUR2_s3,main=paste('OUR2 (s=3) = ',signif(err_s[3],digits = 4),sep=''))
    }
  }
  
  if(method$OUR3){
    ###   use mds to obtain a s-dim version of each data point, use SCMS to estimate the manifold
    ###   and estimation geo with robust isomap.
    ###   bandwidth is selected with Silverman's rule of thumb
    print("mds of smooth data + scms_h_heuri + robust isomap (no oracle)")
    
    err_s<- rep(0,3)
    for(s in 1:3){
      # Projection of the data in dim s using mds
      euc_dis_mat=dist(smooth_data)
      mds_mat = cmdscale(euc_dis_mat,eig=TRUE, k=s)
      projected = mds_mat$points
      denoised = scms$scms(projected)
      estim_geo_mds_scms_RI<-robust_iso(denoised)
      list_to_return[[paste("estim_geo_OUR3_s",s,sep='')]] = estim_geo_mds_scms_RI
      if(Analytic_geo_available) {
        err<- sqrt(sum((estim_geo_mds_scms_RI -true_geo)^2))/norm_true_geo
        err_s[s]<-err
      }
      
    }
    
    if(plot_true){
      
      par(mfrow=c(1,3))
      image.plot(list_to_return$estim_geo_OUR3_s1,main=paste('OUR3 (s=1) = ',signif(err_s[1],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_OUR3_s2,main=paste('OUR3 (s=2) = ',signif(err_s[2],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_OUR3_s3,main=paste('OUR3 (s=3) = ',signif(err_s[3],digits = 4),sep=''))
    }
  }
  
  
  ###  use random projection to obtain a s-d version of each data point, use SCMS to estimate the manifold
  ###    and estimation geo with Floyd's algo.
  if(method$RP){
    print("Random proj. + scms + Floyd's Algorithm")
    scms_h=seq(0.1,2,by=0.2)
    err_RP_s<- rep(0,3)
    for(s in 1:3){
      geo_proj<-matrix(0,nrow=nb_proj,ncol=samplesize*(samplesize-1)/2)
      op_h_proj<-rep(0,nb_proj)
      op_num_nei_proj<-rep(0,nb_proj)
      err_proj<-rep(0,nb_proj)
    
      for(proj in 1:nb_proj){
        # use random projection to represent data in dim s
        randprojmat = matrix(rnorm(s*K,mean=0,sd=1/sqrt(s)), nrow=s, ncol=K)
        projected = randprojmat %*% t(smooth_data)
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
      estim_geo_RP_scms <- ens_geo_mat + t(ens_geo_mat)
      list_to_return[[paste("estim_geo_RP_s",s,sep='')]] = estim_geo_RP_scms
      err_RP_s[s]<-sqrt(sum((estim_geo_RP_scms -true_geo )^2))/norm_true_geo
    }
    if(plot_true){
      par(mfrow=c(1,3))
      image.plot(list_to_return$estim_geo_RP_s1,main=paste('err RP + scms (s=1) = ',signif(err_RP_s[1],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_RP_s2,main=paste('err RP + scms (s=2) = ',signif(err_RP_s[2],digits = 4),sep=''))
      image.plot(list_to_return$estim_geo_RP_s3,main=paste('err RP + scms (s=3) = ',signif(err_RP_s[3],digits = 4),sep=''))
    }
  }
  
  if(method$L2){
    ## L2 distance from smooth data
    
    print("pairwise L2 distances of smooth data")
    
    estim_L2<-as.matrix(dist(smooth_data,diag=TRUE,upper=TRUE))
    
    if(plot_true){
      par(mfrow=c(1,1))
      image.plot(estim_L2,main='L2 distances')
    }

    # Update list of estimation geo
    list_to_return[["estim_L2"]] <- estim_L2
  }
  
  if(method$w_L2){
    ## weigthed L2 distance from smooth data
    
    temp<- pairwise_weigthed_L2(discrete_data,grid[1,],reg_grid)
    
    # nbas=15
    # nbasis.w=10
    # 
    # print("pairwise weigthed L2 distances of smooth data")
    # y<- t(discrete_data)/(sqrt((b-a)/K))
    # bsb = create.bspline.basis(range(reg_grid), nbasis=nbas)
    # B = eval.basis(reg_grid, bsb)
    # P = getbasispenalty(bsb)
    # delta=1/K
    # 
    # arrayVp=array(NA, c(nbas,nbas, samplesize))
    # coef=matrix(NA, nbas, samplesize)
    # 
    # for (i in 1:samplesize){
    #   swmod = gam(y[,i]~B-1,paraPen=list(B=list(P)), method="REML")
    #   arrayVp[,,i]=swmod$Vp
    #   coef[,i]=swmod$coefficients
    # }
    # 
    # 
    # bsb.w = create.bspline.basis(range(reg_grid), nbasis=nbasis.w)
    # B.grid.w=eval.basis(reg_grid, bsb.w)
    # 
    # fit=weight.minCV(coef=coef, arrayVp=arrayVp, B.grid=B, B.grid.weight=B.grid.w, t.grid=reg_grid)
    # weight=fit$weight
    # #plot(reg_grid, weight, type='l')
    estim_w_L2<-temp$pairwise_L2
    
    if(plot_true){
      par(mfrow=c(1,1))
      image.plot(estim_w_L2,main='weigthed L2 distances')
    }
    
    # Update list of estimation geo
    list_to_return[["estim_w_L2"]] <- estim_w_L2
  }
  
  return(list_to_return)

}

##### Smooth data one curve at the time with bspline and penatlty on second derivative

smooth_FD_bspline <- function(discrete_data,grid,reg_grid){
  
  K<- ncol(discrete_data)
  nbasis <- K+2
  basis_obj <- create.bspline.basis(c(reg_grid[1],reg_grid[K]),nbasis)
  loglambda <- seq(-10,-1,by=0.5)
  gcv<-rep(0,length(loglambda))
  for(i in 1:length(loglambda)){
    fd_par_obj <- fdPar(fdobj=basis_obj,Lfdobj=2,lambda=10^loglambda[i])
    smooth_res <- smooth.basis(t(grid),t(discrete_data),fd_par_obj)
    gcv[i]=mean(smooth_res$gcv)
  }
  ind_m=which(gcv==min(gcv))
  fd_par_obj <- fdPar(fdobj=basis_obj,Lfdobj=2,lambda=10^loglambda[ind_m])
  smooth_res <- smooth.basis(t(grid),t(discrete_data),fd_par_obj)
  smooth_curve_fd <- smooth_res$fd
  #smooth_curve <- eval.fd(reg_grid,smooth_curve_fd)
  
  return(smooth_curve_fd)
  
  #return(t(smooth_curve))
  
}