##### Dataset generators

# Scenario 1 : random function from the manifold M={X(t)=u(h_alpha(t)},with h_alpha(t)=t-alpha and u(t)=t
# TODO Marie: the analytic geodesic distance calculation is incorrect for scenario 2
# Scenario 2 : random function from the manifold M={X(t)=density N(alpha,1)}

# TODO Marie/Susan: we need more functional scenarios!

## Input 
# samplesize
# K = number of grid points
# SNR (signal to noise ratio)
# reg_sampling : if 0 points on the manifold are draw from a norm dist. and if 1 regular sampling 
# a,b : functions X(t) are defined for t in [a,b]
# com_grid : if 1 each curve is observed on the same grid and if 0 the grid for each curve is randomly genarated from unif[a,b]
# plot_true : if 1 plot the true data, the observed data and the true geodesic matrix

## Output
# noiseless_data : samplesize x K matrix containing the original data (no noise)
# noisy_data : samplesize x K matrix containing the observed data
# analytic_geo : samplesize x samplesize matrix containing analytic pairwise geo disctance
# grid : samplesize x K matrix containing the grid on which each data is observed
# reg_grid : vector of dim K containing a common grid to use for smoothing

#Example of call
# data<- sim_functional_data(1,20,15,-3,3,0.5,1,1,1)

library(fields)
source('full_geo_from_adj_geo.R')

sim_functional_data<-function(sce,samplesize=100,K=20,a=0,b=1,SNR=1,reg_sampling=1,com_grid=1,plot_true=1){
  
  if(sce == 1){
    if(reg_sampling==0){
      Z <- rnorm(samplesize,0,0.3)
      alpha <- apply(cbind(rep(-1,samplesize),Z),1,max)
      alpha=sort(alpha)
    } else if(reg_sampling==1){
      alpha <- seq(-0.9,3,length.out=samplesize)
    }
    
    mu_t <- function(t,al){
      h_t <- t-al
    }
    
    adja_geo <- (alpha[-1]- alpha[-samplesize])*sqrt(b-a)
    ### Calculate the analytic geodesic matrix
    analytic_geo <- full_geo(adja_geo,samplesize)
    
  } else if(sce ==2){
    
    if(reg_sampling==0){
      Z <- rnorm(samplesize,0,0.3)
      alpha <- apply(cbind(rep(-1,samplesize),Z),1,max)
      alpha=sort(alpha)
    } else if(reg_sampling==1){
      alpha <- seq(-2,2,length.out=samplesize)
    }
    
    mu_t <- function(t,al){
      fct <- dnorm(t,al,1)
    }
    
    adja_geo <- (alpha[-1]- alpha[-samplesize])/(2*pi^(1/4))
    ### Calculate the analytic geodesic matrix
    analytic_geo <- full_geo(adja_geo,samplesize)
    
  }
  
  noiseless_data <- matrix(ncol=K,nrow=samplesize)
  grid<-matrix(ncol=K,nrow=samplesize)
  reg_grid=seq(a,b,length.out=K)
  if(com_grid==0){
    reg_noiseless_data = matrix(ncol=K,nrow=samplesize)
    for(i in 1:samplesize){
      tmp_grid=sort(runif(K,a,b))
      grid[i,]=tmp_grid
      noiseless_data[i,] <- mu_t(tmp_grid,alpha[i])
      reg_noiseless_data[i,]<-mu_t(reg_grid,alpha[i])
    }
  }else if (com_grid==1){
    for(i in 1:samplesize){
      grid[i,]=reg_grid
      noiseless_data[i,] <- mu_t(reg_grid,alpha[i])
    }
    reg_noiseless_data=noiseless_data
  }
    
  mean_signal= apply(reg_noiseless_data,2,mean)
  var_signal= (1/(samplesize*K))*sum((reg_noiseless_data-matrix(mean_signal,ncol=K,nrow=samplesize,byrow=TRUE))^2)
  sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
  epsilon<- matrix(rnorm(samplesize*K,0,sd_noise),nrow=samplesize)
  noisy_data <- (noiseless_data + epsilon)
    
  if (plot_true==1){
    par(mfrow=c(1,3))
    matplot(t(grid),t(noiseless_data),main="True data",type='l', col=rainbow(samplesize))
    matplot(t(grid),t(noisy_data),main="Observed data",type='l', col=rainbow(samplesize))
    image.plot(analytic_geo,main='analytic geodesic')
  }
    
  return(list('noiseless_data'=noiseless_data,'noisy_data'=noisy_data,'analytic_geo'=analytic_geo,'grid'=grid,'reg_grid'=reg_grid))
}