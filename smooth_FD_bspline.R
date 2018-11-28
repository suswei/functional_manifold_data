##### Smooth data one curve at the time with bspline and penatlty on second derivative

smooth_FD_bspline <- function(discrete_data,grid,reg_grid){
  
  K = ncol(discrete_data)
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
  smooth_curve_tmp <- smooth_res$fd
  smooth_curve <- eval.fd(reg_grid,smooth_curve_tmp)
  return(t(smooth_curve))
  
}