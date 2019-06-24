#' ---
#' title: Functional manifold data examples
#' ---

#' \newcommand {\To}{\rightarrow}
#'\newcommand {\TO}{\Rightarrow}
#'\newcommand {\R}{\mathbb{R}}
#'\newcommand {\Prob}{\mathbb{P}}
#'\newcommand{\E}{\mathbb{E}}
#'\newcommand {\cov}{\textrm{Cov}}
#'\newcommand {\var}{\textrm{Var}}
#'\newcommand {\1}{\textrm{\textbf{1}}}
#'\newcommand{\M}{\mathcal{M}}
#'
#' ##Scenario 1: 
#' Consider the manifold $$ \M = \{X_\alpha: -1 \le \alpha \le 1\} $$ with the $L_2$ inner product as the metric tensor,
#' where $X_\alpha: [a,b] \to \R$ is given by $X_\alpha(t) = t-\alpha$. Below, we set $a=-4$ and $b=4$.
#' The shortest path $\gamma: [0,1] \to \M$ between two functions $X_{\alpha_1} = (t-\alpha_1)$ and $X_{\alpha_2} = (t-\alpha_2)$ in $\M$ is clearly given by $\gamma(t) = t X_{\alpha_1} + (1-t) X_{\alpha_2}$.
#' The geodesic distance between $X_{\alpha_1}$ and $X_{\alpha_2}$ is $$L(\gamma) = \int_0^1 || \dot \gamma(t) ||_{L^2} \,dt = ||X_{\alpha_1} - X_{\alpha_2}||_{L_2} = (\alpha_1-\alpha_2)\sqrt{b-a}.$$ 
#'
#'
#' ##Scenario 2: 
#' This scenario is modified from what is referred to as Manifold 2 in Chen and Muller 2012 by fixing the variance of the normal density to be $1$. 
#' We have $$\M =  \left \{X_\beta: \beta \in [-1,1], t \in [a,b]\right \}$$ with the $L_2$ inner product as the metric tensor of $\M$, 
#' where $X_\beta: [a,b] \to \R$ is given by $X_\beta(t) = \frac{1}{\sqrt{2\pi}} \exp{[-\frac{1}{2}(t-\beta)^2]}$. Below, we set $a=-4$ and $b=4$.
#' The geodesic distance between the curves $X_{\beta_1}$ and $X_{\beta_2}$ is given by
#'\begin{eqnarray*}
#'d(X_{\beta_1},X_{\beta_2}) &=& \int_{\beta_1}^{\beta_2} \left \| \frac{d X_\beta (t)}{d\beta} \right\|_{L^2} d\beta \\
#'&=&  \int_{\beta_1}^{\beta_2} \sqrt{ \frac{1}{2\sqrt{\pi}} \int_{-4}^4 \frac{1}{\sqrt{\pi}} \exp\{-(t-\beta)^2\}(t-\beta)^2 dt  }  \  d\beta \\
#'&=&   \int_{\beta_1}^{\beta_2} \sqrt{ \frac{1}{2\sqrt{\pi}} \int_{-4}^4(t-\beta)^2 f(t) dt  }  \  d\beta, \textrm{ where $f$ is the density of a N$(\beta,1/2)$ }   \\
#'&\approx&  \int_{\beta_1}^{\beta_2}  \sqrt{ \frac{1}{2\sqrt{\pi}} \frac{1}{2}} \ d\beta \\
#'&=& (\beta_2-\beta_1) \frac{1}{2\pi^{1/4}},
#'\end{eqnarray*}
#' where the approximation comes from the fact that we are integrating on $[a,b]=[-4,4]$ and not on $\R$. 
#' We can see this manifold is isometric, since the geodesic distance between $X_{\beta_1}$ and $X_{\beta_2}$ in $\M$ is the Euclidan distance between the $\beta$'s, up to some scaling factor. 
#' Note that the "straight" line connecting $X_{\beta_1}$ and $X_{\beta_2}$ in $\M$ does not always stay inside of $\M$, so we cannot employ the calculation technique of Scenario 1. 
#'
#' ##Scenario 3: 
#' Fix $\mu_1,\sigma_1^2, \mu_2, \sigma_2^2$. Let $f_1$ be the normal density with mean $\mu_1$ and variance $\sigma_1^2$. Define $f_2$ analogously.
#' Consider the manifold $$ \M = \{ X_c: 0 \le c \le 1 \} $$ with the $L_2$ inner product as the metric tensor where
#' $X_c: [a,b] \to \R$ is given by $X_c(t) = c f_1(t) + (1-c) f_2(t)$.
#' Again, since $$ \gamma(t) := t X_{c_1} + (1-t) X_{c_2}$$ belongs to $M$ for all $t$, it must be the shortest path between $X_1$ and $X_2$. 
#' The geodesic distance between $X_{c_1}$ and $X_{c_1}$ is
#' $$\int_0^1 || \dot \gamma(t) ||_{L^2} \,dt = ||X_1 - X_2||_{L_2} = (c_1-c_2) ||f_1 - f_2||_{L_2}$$ which has no easy analytic solution but which we estimate using numerical integration. 
#' One may also verify the equation 
#' $$ \int_{c_1}^{c_2} \left \| \frac{d X_c (t)}{dc} \right\|_{L^2} dc $$
#' gives the same result.
#'
#'
#' ##Scenario 4: 
#' It was shown in [Srivastava2007](https://www-sop.inria.fr/ariana/Projets/Shapes/ThirdYearReport/JoshietalCVPR07b.pdf) that the square root representation of probability density functions has a nice closed form geodesic.
#' They consider the manifold $$ \M = \{ \psi:[0,1] \to \R : \psi \ge 0, \int_0^1 \psi^2(s) \,ds = 1 \}$$ with the metric tensor given by the Fisher-Rao metric tensor
#' $$ <v_1,v_2> = \int_0^1 v_1(s) v_2(s) \,ds $$ for two tangent vectors $v_1,v_2 \in T_\psi(\M)$.
#' Note that this concides with the $L_2[0,1]$ inner product.
#' [Srivastava2007] showed that the geodesic distance between any two $\psi_1$ and $\psi_2$ in $\M$ is simply
#' $$d(\psi_1,\psi_2) = \cos^{-1}<\psi_1,\psi_2>$$.
#' We will specifically examine the square root of $Beta(\alpha,\beta)$ distributions which is supported on $[0,1]$. That is, 
#' $$ M = \{ \psi_{\alpha,\beta}: 1 \le \alpha \le 5, 2 \le \beta \le 5\} $$
#' where $\psi_{\alpha,\beta}: [0,1] \to \R$ is the pdf of $Beta(\alpha,\beta)$.


## Input 
# sce: scenario number
# samplesize: number of functional data 
# K: number of time grid points
# SNR (signal to noise ratio)
# reg_sampling : 0 if manifold parameters are drawn randomly and 1 if deterministically
# com_grid : if 1 each curve is observed on the same grid and if 0 the grid for each curve is randomly genarated from unif[a,b]
# plot_true : if 1 plot the true data, the observed data and the true geodesic matrix

## Output
# noiseless_data : samplesize x K matrix containing the original data (no noise)
# noisy_data : samplesize x K matrix containing the observed data
# analytic_geo : samplesize x samplesize matrix containing analytic pairwise geo disctance
# grid : samplesize x K matrix containing the grid on which each data is observed
# reg_grid : vector of dim K containing a common grid to use for smoothing



#Example of call
# data<- sim_functional_data(1,100,15,0.5,1,1,1)

library(fields)
source('full_geo_from_adj_geo.R')

sim_functional_data<-function(sce,samplesize=100,K=30,SNR=1,reg_sampling=1,com_grid=1,plot_true=1){
  if(sce == 1){
    a<- -4
    b<- 4
    if(reg_sampling==0){
      alpha<- runif(samplesize,-1,1)
      alpha=sort(alpha)
    } else if(reg_sampling==1){
      alpha <- seq(-1,1,length.out=samplesize)
    }
    
    mu_t <- function(t,al){
      h_t <- t-al
    }
    
    # TODO: Add mathematical derviation for adja_geo
    adja_geo <- (alpha[-1]- alpha[-samplesize])*sqrt(b-a)
    ### Calculate the analytic geodesic matrix
    analytic_geo <- full_geo(adja_geo,samplesize)
    
  } else if(sce ==2){
    a<- -4
    b<- 4
    if(reg_sampling==0){
      alpha<- runif(samplesize,-1,1)
      alpha=sort(alpha)
    } else if(reg_sampling==1){
      alpha <- seq(-1,1,length.out=samplesize)
    }
    
    mu_t <- function(t,al){
      fct <- dnorm(t,al,1)
    }
    
    adja_geo <- (alpha[-1]- alpha[-samplesize])/(2*pi^(1/4))
    ### Calculate the analytic geodesic matrix
    analytic_geo <- full_geo(adja_geo,samplesize)
  }else if(sce==3){
    a <- -4
    b <- 4
    
    # we make samplesize different combinations of the parameters alpha and sigma
    nb_alpha<- 10
    nb_beta<- samplesize/10
    if(reg_sampling==0){
      alpha<- runif(samplesize,-1,1)
      alpha=sort(alpha)
      sig<- runif(nb_beta,0.5,1.5)
      sig<-sort(sig)
    } else if(reg_sampling==1){
      alpha <- seq(-1,1,length.out=nb_alpha)
      sig<- seq(0.5,1.5,length.out=nb_beta)
    }
    alpha_beta<- expand.grid(alpha,sig)
    
    mu_t <- function(t,alpha_beta){
      fct <- dnorm(t,alpha_beta[,1],alpha_beta[,2])
    }
    
    analytic_geo <- matrix(0,samplesize,samplesize)
    grid_int <- seq(a,b,0.01)
    for(comb1 in 1:(samplesize-1)){
      for(comb2 in comb1:(samplesize)){
        
        f <- function(x){
          (dnorm(x,alpha_beta[comb1,1],alpha_beta[comb1,2]) - dnorm(x,alpha_beta[comb2,1],alpha_beta[comb2,2]))^2
        }
        
        analytic_geo[comb1,comb2]<-  sqrt(integrand)
      }
    }
    analytic_geo <- analytic_geo + t(analytic_geo)
    
  }else if(sce==4){
    
    a = 0
    b = 1
    
    nb_alpha<- 10
    nb_beta<- samplesize/10
    if(reg_sampling==0){
      alpha<- runif(samplesize,1,5)
      alpha=sort(alpha)
      beta<- runif(nb_beta,2,5)
      beta<-sort(sig)
    } else if(reg_sampling==1){
      alpha <- seq(1,5,length.out=nb_alpha)
      beta<- seq(2,5,length.out=nb_beta)
    }
    alpha_beta<- expand.grid(alpha,beta)
    
    mu_t <- function(t,alpha_beta){
      fct <- sqrt(dbeta(t,alpha_beta[,1],alpha_beta[,2]))
    }
    
    analytic_geo <- matrix(0,samplesize,samplesize)
    grid_int <- seq(a,b,0.01)
    for(comb1 in 1:(samplesize-1)){
      for(comb2 in (comb1+1):samplesize){
        
        f <- function(x){
          sqrt(dbeta(x,alpha_beta[comb1,1],alpha_beta[comb1,2])) * sqrt(dbeta(x,alpha_beta[comb2,1],alpha_beta[comb2,2]))
        }
        
        analytic_geo[comb1,comb2]<-  
          acos(integrate(f,lower=a,upper=b)$value)
      }
    }
    analytic_geo <- analytic_geo + t(analytic_geo)
  }
  
  noiseless_data <- matrix(ncol=K,nrow=samplesize)
  grid<-matrix(ncol=K,nrow=samplesize)
  reg_grid=seq(a,b,length.out=K)
  if(com_grid==0){
    reg_noiseless_data = matrix(ncol=K,nrow=samplesize)
    for(i in 1:samplesize){
      tmp_grid=sort(runif(K,a,b))
      grid[i,]=tmp_grid
      if(sce==1 || sce==2){
        noiseless_data[i,] <- mu_t(tmp_grid,alpha[i])
        reg_noiseless_data[i,]<-mu_t(reg_grid,alpha[i])
      } else if(sce==3 || sce==4){
        noiseless_data[i,] <- mu_t(tmp_grid,alpha_beta[i,])
        reg_noiseless_data[i,]<-mu_t(reg_grid,alpha_beta[i,])
      }
      
    }
  }else if (com_grid==1){
    for(i in 1:samplesize){
      grid[i,]=reg_grid
      if(sce==1 || sce==2){
        noiseless_data[i,] <- mu_t(reg_grid,alpha[i])
      }else if(sce==3 || sce==4){
        noiseless_data[i,] <- mu_t(reg_grid,alpha_beta[i,])
      }
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

