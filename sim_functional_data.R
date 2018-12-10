##### Dataset generators

# Scenario 1 : random function from the manifolg M={X(t)=u(h_alpha(t)},with h_alpha(t)=t-alpha and u(t)=t
# TODO Marie: the true geodesic distance calculation is incorrect for scenario 2
# Scenario 2 : random function from the manifolg M={X(t)=density N(alpha,1)}

# TODO Marie/Susan: we need more functional scenarios!

## Input 
# samplesize
# K = number of grid points
# SNR (signal to noise ratio)
# reg_sampling : if 0 points on the manifold are draw from a norm dist. and if 1 regular sampling 
# a,b : functions X(t) are defined for t in [a,b]
# com_grid : if 1 each curve is observed on the same grid and if 0 the grid for each curve is randomly genarated from unif[a,b]
# plot_true : if 1 plot the true data, the observed the data and the true geodesic matrix

## Output
# true_data : samplesize x K matrix containing the original data (no noise)
# discrete_data : samplesize x K matrix containing the observed data
# true_geo : samplesize x samplesize matrix containing true pairwise geo disctance
# grid : samplesize x K matrix containing the grid on which each data is observed
# reg_grid : vector of dim K containing a common grid to use for smoothing

#Example of call
# data<- sim_functional_data(1,20,15,-3,3,0.5,1,1,1)

sim_functional_data<-function(sce,samplesize,K,a,b,SNR,reg_sampling,com_grid,plot_true){
  
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
    
    ### Calculate the true geodesic matrix
    adja_geo <- (alpha[-1]- alpha[-samplesize])*sqrt(b-a)
    true_geo <- full_geo(adja_geo,samplesize)
    
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
    
    ### Calculate the true geodesic matrix
    adja_geo <- (alpha[-1]- alpha[-samplesize])/(2*pi^(1/4))
    true_geo <- full_geo(adja_geo,samplesize)
    
  }
  
  true_data <- matrix(ncol=K,nrow=samplesize)
  grid<-matrix(ncol=K,nrow=samplesize)
  reg_grid=seq(a,b,length.out=K)
  if(com_grid==0){
    reg_true_data = matrix(ncol=K,nrow=samplesize)
    for(i in 1:samplesize){
      tmp_grid=sort(runif(K,a,b))
      grid[i,]=tmp_grid
      true_data[i,] <- mu_t(tmp_grid,alpha[i])
      reg_true_data[i,]<-mu_t(reg_grid,alpha[i])
    }
  }else if (com_grid==1){
    for(i in 1:samplesize){
      grid[i,]=reg_grid
      true_data[i,] <- mu_t(reg_grid,alpha[i])
    }
    reg_true_data=true_data
  }
    
  mean_signal= apply(reg_true_data,2,mean)
  var_signal= (1/(samplesize*K))*sum((reg_true_data-matrix(mean_signal,ncol=K,nrow=samplesize,byrow=TRUE))^2)
  sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
  epsilon<- matrix(rnorm(samplesize*K,0,sd_noise),nrow=samplesize)
  discrete_data <- (true_data + epsilon)
    
  if (plot_true==1){
    par(mfrow=c(1,3))
    matplot(t(grid),t(true_data),main="True data",type='l', col=rainbow(samplesize))
    matplot(t(grid),t(discrete_data),main="Observed data",type='l', col=rainbow(samplesize))
    image.plot(true_geo,main='True geodesic')
  }
    
  return(list('true_data'=true_data,'discrete_data'=discrete_data,'true_geo'=true_geo,'grid'=grid,'reg_grid'=reg_grid))
}