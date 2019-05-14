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
#' Consider the manifold $$ \M = \{X(t) = t-\alpha, -0.9 \le \alpha \le 3, t \in [a,b]\} $$ with the $L_2$ inner product as the metric tensor.
#' The shortest path $\gamma: [0,1] \to \M$ between two functions $X_1 = (t-\alpha_1)$ and $X_2 = (t-\alpha_2)$ in $\M$ is clearly given by $\gamma(t) = t X_1 + (1-t) X_2$.
#' The geodesic distance between $X_1$ and $X_2$ is $$L(\gamma) = \int_0^1 || \dot \gamma(t) ||_{L^2} \,dt = ||X_1 - X_2||_{L_2} = (\alpha_1-\alpha_2)\sqrt{b-a}.$$ 
#'
#'
#' ##Scenario 2: 
#' This is based on Manifold 2 in Chen and Muller 2012 where we fix the variance to be $1$. The paper claims that in this case the resulting manifold $$\M =  \left \{ X_\beta \in L^2([-4,4]) : X_\beta(t) = \frac{1}{\sqrt{2\pi}} \exp{[-\frac{1}{2}(t-\beta)^2]}, \beta \in \R\right \}.$$ is isometric.
#' Chen and Muller are also working with the $L_2$ inner product as their metric tensor of $\M$. 
#' Isometric means the geodesic distance between $X_1$ and $X_2$ in $\M$ is the Euclidan distance between the $\alpha$'s up to some scaling factor. ???Why is there a scaling factor???
#' Note that the "straight" line connecting $X_1$ and $X_2$ in $\M$ does not always stay inside of $\M$ since $t N(\alpha_1,1) + (1-t) N(\alpha_2,1)$ no longer has variance $1$ for all $t \in [0,1]$. 
#'  The geodesic distance between the curves $X_{\beta_1}$ and $X_{\beta_2}$ is 
#'\begin{eqnarray*}
#'d(X_{\beta_1},X_{\beta_2}) &=& \int_{\beta_1}^{\beta_2} \left \| \frac{d X_\beta (t)}{d\beta} \right\|_{L^2} d\beta \\
#'&=&  \int_{\beta_1}^{\beta_2} \sqrt{ \frac{1}{2\sqrt{\pi}} \int_{-4}^4 \frac{1}{\sqrt{\pi}} \exp\{-(t-\beta)^2\}(t-\beta)^2 dt  }  \  d\beta \\
#'&=&   \int_{\beta_1}^{\beta_2} \sqrt{ \frac{1}{2\sqrt{\pi}} \int_{-4}^4(t-\beta)^2 f(t) dt  }  \  d\beta, \textrm{ where $f$ is the density of a N$(\beta,1/2)$ }   \\
#'&\approx&  \int_{\beta_1}^{\beta_2}  \sqrt{ \frac{1}{2\sqrt{\pi}} \frac{1}{2}} \ d\beta \\
#'&=& (\beta_2-\beta_1) \frac{1}{2\pi^{1/4}},
#'\end{eqnarray*}
#'where the approximation comes from the fact that we are integrating on $[-4,4]$ and not on $\R$. Moreover the mean $\beta$ has to be close to 0 for the approximation to make sense.
#'
#' #Scenario 3: ???This needs to be implemented eventually???
#' Consider the manifold $$ \M = \{ \text{probability density function } N(\alpha,\sigma^2): \alpha \in \mathbb R, \sigma > 0 \} $$ with the $L_2$ inner product as the metric tensor.
#' Let $X_1$ be the pdf of $N(\alpha_1,\sigma_1^2)$ and $X_2$ be the pdf of $N(\alpha_2,\sigma_2^2)$. 
#' Again, since $$ \gamma(t) := t X_1 + (1-t) X_2$$ belongs to $M$ for all $t$, it must be the shortest path between $X_1$ and $X_2$. The geodesic distance between $X_1$ and $X_2$ is
#' $$\int_0^1 || \dot \gamma(t) ||_{L^2} \,dt = ||X_1 - X_2||_{L_2}$$ which has no easy analytic solution.
#'
#'
#' #Scenario 4: ???This needs to be implemented eventually???
#' This is taken from [Srivastava2007](https://www-sop.inria.fr/ariana/Projets/Shapes/ThirdYearReport/JoshietalCVPR07b.pdf). 
#' Consider the manifold $$ \M = \{ \psi:[0,1] \to \R : \psi \ge 0, \int_0^1 \psi^2(s) \,ds = 1 \}$$ with the Fisher-Rao metric tensor given by
#' $$ <v_1,v_2> = \int_0^1 v_1(s) v_2(s) \,ds $$ for two tangent vectors $v_1,v_2 \in T_\psi(\M)$. The geodesic distance between any two $\psi_1$ and $\psi_2$ in $\M$ is simply
#' $$d(\psi_1,\psi_2) = \cos^{-1}<\psi_1,\psi_2>$$.

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

sim_functional_data<-function(sce,samplesize=100,K=30,a=-4,b=4,SNR=1,reg_sampling=1,com_grid=1,plot_true=1){
  
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
      alpha <- seq(-1,1,length.out=samplesize)
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

  
