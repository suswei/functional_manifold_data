# generate 2D manifold and noisy observations off it
# Input
#   samplesize: self explanatory
#   SNR: signal to noise ratio SNR= 10*log_10 (var_signal/sigma^2), where sigma^2 = variance of noise (ref. Orientability and Diffusion Maps, Amit Singera and Hau-tieng Wu)
#   plotTrue: plot of real and noisy manifold
#   reg_sampling: TRUE = regular and FALSE = random uniformly
#   min_ana_num: if calculus of true geodesic done numerically, = number of points used in approx
#   name options:
#     archimedean-spiral: complete
#     sin-curve: 
#     cross: complete
#     circle: complete
#     right-angle: complete
#     half-moon (circle of radius 1 with unregular sampling) :complete
#     bananas (two unconnected components) : complete (geo calculated by numeric integration)
#     line: complete
#     sin-cos-curve: complete (geo calculated by numeric integration)
#     hetero-spiral: complete (geo calculated by numeric integration)

# OUTPUT
# data: samplesize x 2 matrix of noisy manifold observations
# true_mani: samplesize x 2 matrix of true manifold observations
# true_gep: samplesize x samplesize matrix of true geodesic distances
# plots true and noisy manifold observations


EuclideanExamples <- function(name, samplesize, SNR, plotTrue,reg_sampling,min_ana_num){
  
  if ( name == "archimedean-spiral") {
    
    if(reg_sampling==TRUE){
      t = seq(0, 5*pi, length.out = samplesize)
    }else {
      t = runif(samplesize,0, 5*pi)
      t = sort(t)
    }
    
    x = 1 * t * cos(t)
    y = 1 * t * sin(t)
    true_mani = cbind(x, y)
    
    adja_geo <- arc_length_spiral(t[-1])- arc_length_spiral(t[-samplesize])
    true_geo <- full_geo(adja_geo,samplesize)
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(length(t), sd = sd_noise)
    e2 = rnorm(length(t), sd = sd_noise)
    data = cbind( x+e1 , y+e2 )
    
    
  } else if ( name == "sin-curve") {
    
    # calling manifold_data.py will return double the samplesize,
    # first set is the noisy manifold data, second set is the true manifold
    source_python('manifold_data.py')
    stacked = manifold_data(samplesize)
    
    true_mani = stacked[-(1:samplesize), ]
    data = stacked[1:samplesize, ]
    
    
  } else if ( name == "cross") {
    
    if(reg_sampling==TRUE){
      x <- c(seq(-4,4,length.out=samplesize/2),rep(0,samplesize/2))
      y <- c(rep(0,samplesize/2),seq(-4,4,length.out=samplesize/2))
    }else {
      x <- c(sort(runif(samplesize/2,-4,4)),rep(0,samplesize/2))
      y <- c(rep(0,samplesize/2),sort(runif(samplesize/2,-4,4)))
    }
    
    true_mani <- cbind(x,y)
    
    true_geo <- matrix(0,ncol=samplesize,nrow=samplesize)
    for(i in 1:(samplesize-1)){
      for(j in (i+1):samplesize){
        true_geo[i,j] <- abs(true_mani[i,1]-true_mani[j,1])+abs(true_mani[i,2]-true_mani[j,2])
      }
    }
    true_geo <- true_geo+t(true_geo)
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(samplesize,sd=sd_noise)
    e2 = rnorm(samplesize, sd=sd_noise)
    data = cbind(x+e1,y+e2)
    
    
    
    
  } else if ( name == "circle") {
    
    if(reg_sampling==TRUE){
      angle = seq( 0,2*pi,length.out=samplesize)
    }else {
      angle = runif(samplesize, 0,2*pi)
      angle <- sort(angle)
    }
    
    radius <- 10
    x=radius*cos(angle)
    y=radius*sin(angle)
    true_mani = cbind(x ,y)

    
    true_geo <- matrix(0,ncol=samplesize,nrow=samplesize)
    for(i in 1:(samplesize-1)){
      for(j in (i+1):samplesize){
        true_geo[i,j] <- min(abs(angle[i]-angle[j]),2*pi-abs(angle[i]-angle[j]))*radius		
      }
    }
    true_geo <- true_geo+t(true_geo)
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(samplesize,sd=sd_noise)
    e2 = rnorm(samplesize, sd=sd_noise)
    data = cbind(x+e1,y+e2)
    
    
  } else if ( name == "right-angle") {
    
    
    if(reg_sampling==TRUE){
      x = c( seq(-10, 0, length.out = samplesize/2), rep(0, samplesize/2) )
      y = c( rep(0, samplesize/2), seq(0, 10, length.out = samplesize/2) )
      
    }else {
      x <- c(sort(runif(samplesize/2,-10,0)),rep(0,samplesize/2))
      y <- c(rep(0,samplesize/2),sort(runif(samplesize/2,0,10)))
    }
    
    true_mani = cbind(x, y)
    
    true_geo <- matrix(0,ncol=samplesize,nrow=samplesize)
    for(i in 1:(samplesize-1)){
      for(j in (i+1):samplesize){
        true_geo[i,j] <- abs(true_mani[i,1]-true_mani[j,1])+abs(true_mani[i,2]-true_mani[j,2])
      }
    }
    true_geo <- true_geo+t(true_geo)
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(samplesize, sd = sd_noise)
    e2 = rnorm(samplesize, sd = sd_noise)
    data = cbind( x+e1 , y+e2 )
    
    
  } else if ( name == "half-moon") {
    
    angle = sort(rnorm(samplesize, mean=0, sd=1))
    radius <- 1
    true_mani = cbind( radius*cos(angle), radius*sin(angle) )
    
    true_geo <- matrix(0,ncol=samplesize,nrow=samplesize)
    for(i in 1:(samplesize-1)){
      for(j in (i+1):samplesize){
        true_geo[i,j] <- min(abs(angle[i]-angle[j]),2*pi-abs(angle[i]-angle[j]))*radius		
      }
    }
    true_geo <- true_geo+t(true_geo)
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(samplesize, sd = sd_noise)
    e2 = rnorm(samplesize, sd = sd_noise)
    data = cbind( x+e1 , y+e2 )
    
    
  } else if ( name == "bananas") {
    
    
    if(reg_sampling==TRUE){
      phi = seq(0, 2*pi,length.out=samplesize)
    }else {
      phi = runif(samplesize,0, 2*pi)
      phi=sort(phi)
    }
    
    xtemp = (10)*cos(phi)
    ytemp = (10)*sin(phi)
    
    x = xtemp - 2*(xtemp>0)
    y = ytemp + 10*(xtemp>0)
    
    ind_1st_comp <- which(xtemp>0)
    dim_1st_comp <- length(ind_1st_comp)
    temp <- which((ind_1st_comp[1:(dim_1st_comp-1)]-ind_1st_comp[2:(dim_1st_comp)])< -1)
    ind_1st_comp <- c(ind_1st_comp[(temp+1):dim_1st_comp],ind_1st_comp[1:temp])
    dim_2st_comp <- samplesize-dim_1st_comp
    x1 <- x[ind_1st_comp]
    x2 <- x[-ind_1st_comp]
    y1 <- y[ind_1st_comp]
    y2 <- y[-ind_1st_comp]
    
    true_mani <- cbind(c(x1,x2),c(y1,y2))
    
    adja_geo1 <- sqrt((x1[-1]-x1[-dim_1st_comp])^2+(y1[-1]-y1[-dim_1st_comp])^2)
    true_geo1 <- full_geo(adja_geo1,dim_1st_comp)
    
    adja_geo2 <- sqrt((x2[-1]-x2[-dim_2st_comp])^2+(y2[-1]-y2[-dim_2st_comp])^2)
    true_geo2 <- full_geo(adja_geo2,dim_2st_comp)
    
    true_geo <- matrix(0,samplesize,samplesize)
    true_geo[1:dim_1st_comp,1:dim_1st_comp] <- true_geo1
    true_geo[(dim_1st_comp+1):samplesize,(dim_1st_comp+1):samplesize] <- true_geo2
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(samplesize, sd = 2)
    e2 = rnorm(samplesize, sd = .5)
    data = cbind( x2+e1 , (y2+e2)/3 )
    
  } else if ( name == "line") {
    
    
    if(reg_sampling==TRUE){
      x = seq(from = -3, to = 3, length.out = samplesize)
    }else {
      x = runif(samplesize,-3,  3)
      x=sort(x)
    }
    
    true_mani = cbind(x, rep(0,samplesize))
    
    adja_geo <- abs(x[-1]-x[-samplesize])
    true_geo <- full_geo(adja_geo,samplesize)
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(samplesize, sd=sd_noise)
    e2 = rnorm(samplesize, sd=sd_noise)
    data = cbind( x + e1, e2 )
    
  } else if ( name == "sin-cos-curve") {
    
    if(reg_sampling==TRUE){
      t = seq(-pi/2,pi/2, length.out=samplesize)
    }else {
      t = runif(samplesize,-pi/2,pi/2)
      t=sort(t)
    }
    
    x = 2*sin(t)
    y = cos(t)+cos(2*t)+cos(6*t)
    true_mani = cbind(x, y)
    
    if(samplesize > min_ana_num){
      adja_geo <- sqrt((x[-1]-x[-samplesize])^2+(y[-1]-y[-samplesize])^2)
      true_geo <- full_geo(adja_geo,samplesize)
    }else{
      adja_geo <- c()
      dense_fac <- ceiling(min_ana_num/samplesize)
      for(i in 1:(samplesize-1)){
        int <- seq(t[i],t[i+1],length.out=dense_fac)
        adja_geo <- c(adja_geo,sum(sqrt((2*sin(int[-1])-2*sin(int[-dense_fac]))^2+(cos(int[-1])+cos(2*int[-1])+cos(6*int[-1])-cos(int[-dense_fac])-cos(2*int[-dense_fac])-cos(6*int[-dense_fac]))^2)))
      }
      true_geo <- full_geo(adja_geo,samplesize)
    }
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(length(t), sd = sd_noise)
    e2 = rnorm(length(t), sd = sd_noise)
    data = cbind( x+e1 , y+e2 )
    
  } else if ( name == "hetero-spiral") {
    
    if(reg_sampling==TRUE){
      t = seq(0,5*pi, length.out = samplesize)
    }else {
      t = runif(samplesize,0,5*pi)
      t=sort(t)
    }
    
    x = 1 * t^(3/2) * cos(t)
    y = 1 * t^(3/2) * sin(t)
    true_mani = cbind(x, y)
    
    if(samplesize > min_ana_num){
      adja_geo <- sqrt((x[-1]-x[-samplesize])^2+(y[-1]-y[-samplesize])^2)
      true_geo <- full_geo(adja_geo,samplesize)
    }else{
      adja_geo <- c()
      dense_fac <- ceiling(min_ana_num/samplesize)
      for(i in 1:(samplesize-1)){
        int <- seq(t[i],t[i+1],length.out=dense_fac)
        adja_geo <- c(adja_geo,sum(sqrt((int[-1]^(3/2)*cos(int[-1])-int[-dense_fac]^(3/2)*cos(int[-dense_fac]))^2+(int[-1]^(3/2)*sin(int[-1])-int[-dense_fac]^(3/2)*sin(int[-dense_fac]))^2)))
      }
      true_geo <- full_geo(adja_geo,samplesize)
    }
    
    mean_signal= apply(true_mani,2,mean)
    var_signal= (1/(samplesize*ncol(true_mani)))*sum((true_mani-matrix(mean_signal,ncol=2,nrow=samplesize,byrow=TRUE))^2)
    sd_noise= sqrt(var_signal/(10^(SNR/10)))
    
    e1 = rnorm(length(t), sd=sd.noise)
    e2 = rnorm(length(t), sd=sd.noise)
    data = cbind( x+e1 , y+e2 )
    
  } 
  
  if(plotTrue){
    par(mfrow=c(1,2))
    plot(true_mani, pch=19, xlab='', ylab='', main = paste("true manifold", name, sep = " "))
    plot(data, pch=19, xlab='', ylab='', main = paste("noisy manifold", name, sep = " "))
  }
  
  result = list(data = data, true_mani = true_mani,true_geo=true_geo)
  return(result)
  
}













