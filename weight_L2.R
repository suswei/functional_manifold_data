
### Weigthed L2 distance from Chen et al (Biometrics)

pairwise_weigthed_L2<- function(discrete_data,grid,reg_grid){
  
  nbas=15
  nbasis.w=10
  samplesize<- nrow(discrete_data)
  y<- t(discrete_data)
  bsb = create.bspline.basis(range(grid), nbasis=nbas)
  B = eval.basis(grid, bsb)
  P = getbasispenalty(bsb)
  delta=reg_grid[2]-reg_grid[1]
  B.grid=eval.basis(reg_grid, bsb)
  
  arrayVp=array(NA, c(nbas,nbas, samplesize))
  coef=matrix(NA, nbas, samplesize)
  smooth_data<- matrix(NA, ncol=length(reg_grid),nrow=samplesize)
  
  for (i in 1:samplesize){
    swmod = gam(y[,i]~B-1,paraPen=list(B=list(P)), method="REML")
    arrayVp[,,i]=swmod$Vp
    coef[,i]=swmod$coefficients
    smooth_data[i,]<- predict(swmod,newdata=list(B=B.grid))
  }
  
  bsb.w = create.bspline.basis(range(grid), nbasis=nbasis.w)
  B.grid.w=eval.basis(reg_grid, bsb.w)
  
  fit=weight.minCV(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.w, t.grid=reg_grid)
  weight=fit$weight
  return(list("pairwise_L2"=as.matrix(sqrt(delta)*dist(t(sqrt(weight)*t(smooth_data)),diag=TRUE,upper=TRUE)),"smooth_data"=smooth_data,"weight"=weight))
  
}




##compute the numerator and denomator matrices B and A respectively
computAB<-function(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=rep(1, ngrid), delta=delta){
  
  nbasis=ncol(B.grid.weight)
  n.curve=ncol(coef)
  A=B1=B2=matrix(0, nbasis, nbasis)
  
  for(i in 1:(n.curve-1)){
    for(j in (i+1):n.curve){
      V=arrayVp[,,i]+arrayVp[,,j]
      varf=diag(B.grid%*%V%*%t(B.grid))
      ebeta=coef[,i]-coef[,j]
      fhat=as.vector(B.grid%*%ebeta)
      
      ### calculate the squared mean
      temp=B.grid.weight*((fhat^2+varf)*sqrt(weight)*delta)
      A1=apply(temp, 2, sum)
      A=A+A1%*%t(A1)
      
      ### calculate the variance
      M=t(B.grid)%*%(B.grid*(weight*delta))
      temp2=diag(B.grid%*%V%*%M%*%V%*%t(B.grid))
      B1=B1+t(B.grid.weight)%*%(B.grid.weight*temp2*delta)
      
      tempv=t(B.grid.weight)%*%(B.grid*(fhat*sqrt(weight)*delta))
      B2=B2+tempv%*%V%*%t(tempv)
    }
  }
  
  Bv=2*B1+4*B2
  
  return(list(A=A, Bv=Bv))
}


########################################################################################
#### coef: an nbasisxn spline coefficient matrix 
#### arrayVp: nbasisXnbasisXn array,  Bayesian variance of the spline coefficient
#### B.grid: ngridXnbasis spline design matrix for the smoothed curves
#### B.grid.weight: spline design matrix for the weight function
#### t.grid: the equally spaced grid points 
########################################################################################

weight.minCV=function(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, t.grid=t.grid, niter=50, tol=1e-3){
  
  ngrid=nrow(B.grid.weight) ###number of grid points
  delta=(max(t.grid)-min(t.grid))/ngrid
  fit=computAB(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=rep(1, ngrid),  delta=delta)
  svdB=svd(fit$Bv)
  Bh=svdB$u%*%diag(svdB$d^(-1/2))%*%t(svdB$v) ##B^{-1/2}
  svdBA=svd(Bh%*%fit$A%*%Bh)
  qnew=Bh%*%svdBA$u[,1]
  CVar2=t(qnew)%*%fit$Bv%*%qnew/(t(qnew)%*%fit$A%*%qnew)
  weight=as.vector((B.grid.weight%*%qnew)^2)
  weight=weight/sum(weight*delta)
  
  #    par(mfrow=c(3,3))
  diff=iter=1
  while(diff>tol & iter<niter){
    cat("iteration=", iter,
        "difference=", diff,
        "CVar^2=", CVar2, "\n")
    
    qold=qnew
    fit=computAB(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=weight, delta=delta)
    svdB=svd(fit$Bv)
    Bh=svdB$u%*%diag(svdB$d^(-1/2))%*%t(svdB$v)
    svdBA=svd(Bh%*%fit$A%*%Bh)
    qnew=Bh%*%svdBA$u[,1]; 
    qnew=(qnew+qold)/2
    weight=as.vector((B.grid.weight%*%qnew)^2);  
    weight=weight/sum(weight*delta)
    diff=mean(abs(B.grid.weight%*%(qold-qnew)))
    iter=iter+1
    CVar2=t(qnew)%*%fit$Bv%*%qnew/(t(qnew)%*%fit$A%*%qnew)
    #        plot(t.grid, weight, pch=20, col='red', ylab='Weight', xlab='t')
  }        
  cat("iteration=", iter,
      "difference=", diff,
      "CVar^2=", CVar2, "\n")
  
  return(list(weight=weight, q=qnew, B.grid.weight=B.grid.weight, CVar2=CVar2))
}




#fit=weight.minCV(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid, t.grid=t.grid)
##
###fit=weight.minCV(coef=coef, arrayVp=arrayVp, B.grid=dB.grid, B.grid.weight=B.grid, t.grid=t.grid)
#fit$q
#fit$weight
#fit$CVar2
#dim(fit$B.grid.weight)
