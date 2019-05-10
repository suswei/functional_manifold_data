#' ## This function generates functional data living in a submanifold of a Hilbert space

#' Scenario 1: Considered is the manifold $$M=\{X(t) = t-\alpha, -0.9 \le \alpha \le 3\} \subset L^2([a,b])$$ with inherited metric from the ambient $L^2[a,b]$ space.
#' The shortest path $\gamma: [0,1] \to M$ between two functions $X_1 = (t-\alpha_1)$ and $X_2 = (t-\alpha_2)$ in $M$ is clearly given by $\gamma(t) = t X_1 + (1-t) X_2$.
#' The geodesic distance between $X_1$ and $X_2$ is $$L(\gamma) := \int_0^1 || \dot \gamma(t) ||_{L^2} \,dt = ||X_1 - X_2||_{L_2} = (\alpha_1-\alpha_2)\sqrt{b-a}.$$ 
#' Scenario 2 : Considered is the manifold $$ M = \{ \text{probability density function} N(\alpha,\sigma^2): \alpha \in \mathbb R, \sigma > 0 \} \subset L^2([-\infty,\infty])$$ with inherited metric from the ambient $L^2[-\infty,\infty]$ space. 
#' Let $X_1$ be the pdf of $N(\alpha_1,\sigma_1^2)$ and $X_2$ be the pdf of $N(\alpha_2,\sigma_2^2)$. 
#' Again, since $$ b\gamma(t) := t X_1 + (1-t) X_2$$ belongs to $M$ for all $t$, it must be the shortest path between $X_1$ and $X_2$. The geodesic distance between $X_1$ and $X_2$ is
#' $$\int_0^1 || \dot \gamma(t) ||_{L^2} \,dt = ||X_1 - X_2||_{L_2}$$ which has no easy analytic solution.
#' Marie's previous example where $$ M = \{ \text{probability density function } N(\alpha,1): \alpha \in \mathbb R \} $$ 
#' is actually much harder to work with because the "straight" line connecting $X_1$ and $X_2$ in $M$ does not always stay inside of $M$. 
#' To see this, note that $t N(\alpha_1,1) + (1-t) N(\alpha_2,1)$ no longer has variance $1$ for all $t \in [0,1]$. 

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
    
    # TODO: Add mathematical derviation for adja_geo
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