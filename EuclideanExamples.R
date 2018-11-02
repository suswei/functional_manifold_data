# TODO Marie: add true geodesic distances
# TODO Marie: add the true manifold to all the examples not labeled complete
# TODO future: not sure if all of these noises are "normal" to the manifold

# generate 2D manifold and noisy observations off it
# Input
#   samplesize: self explanatory
#   sd_noise: for certain simulations, can fiddle with this noise term from outside the function
#   name options:
  #   spiral: complete
  #   sin-curve: complete
  #   cross
  #   circle: complete
  #   right-angle: complete
  #   half-moon: complete
  #   normal
  #   bananas: complete
  #   line: complete
  #   two-clouds:
  #   not-dense-spiral: complete
  #   sin-cos-curve: complete
  #   hetero-spiral: complete

# OUTPUT
# data: samplesize x 2 matrix of noisy manifold observations
# true_mani: samplesize x 2 matrix of true manifold observations
# plots true and noisy manifold observations

library(truncnorm)

EuclideanExamples <- function(name, samplesize, sd_noise, plotTrue){

  if ( name == "spiral") {

    # regular sampling of manifold
    t = seq(pi,5*pi, length.out = samplesize)
    # uniform sampling
    t = runif(samplesize, pi, 5*pi)
    # irregular sampling of manifold
    # t = rtruncnorm(samplesize, a=pi, b=5*pi, mean = 3*pi, sd = 0.8)
    x = 1 * t * cos(t)
    y = 1 * t * sin(t)
    true_mani = cbind(x, y)

    # TODO: is this noise normal to the manifold?
    e1 = rnorm(samplesize, sd=sd_noise)
    e2 = rnorm(samplesize, sd=sd_noise)
    data = cbind(x+e1, y+e2)

  } else if ( name == "sin-curve") {

    # calling manifold_data.py will return double the samplesize,
    # first set is the noisy manifold data, second set is the true manifold
    source_python('manifold_data.py')
    stacked = manifold_data(samplesize)

    true_mani = stacked[-(1:samplesize), ]
    data = stacked[1:samplesize, ]


  } else if ( name == "cross") {

    x = rnorm(samplesize/2, sd = 2)
    y = rnorm(samplesize/2, sd = .2)
    data = data.frame(x=x,y=y)

    x = rnorm(samplesize/2, sd = .2)
    y = rnorm(samplesize/2, sd = 2)
    data = rbind(data, data.frame(x=x,y=y))


  } else if ( name == "circle") {

    angle = runif(samplesize, 0,2*pi)
    true_mani = cbind( 10*cos(angle), 10*sin(angle) )

    # add noise normal to manifold
    alpha = rnorm(samplesize, sd=sd_noise)
    data = cbind( (10+alpha)*cos(angle) , (10+alpha)*sin(angle) )

  } else if ( name == "right-angle") {

    t = seq(0,5*pi, length.out = samplesize)

    x = c( seq(-10, 0, length.out = samplesize/2), rep(0, samplesize/2) )
    y = c( rep(0, samplesize/2), seq(0, 10, length.out = samplesize/2) )
    true_mani = cbind(x, y)

    e1 = rnorm(samplesize, sd = sd_noise)
    e2 = rnorm(samplesize, sd = sd_noise)
    data = cbind( x+e1 , y+e2 )


  } else if ( name == "half-moon") {

    theta = rnorm(samplesize, mean=0, sd=1)
    theta = sort(theta)
    x = cos(theta)
    y = sin(theta)
    true_mani = cbind(x, y)

    e1 = rnorm(samplesize, sd = sd_noise)
    e2 = rnorm(samplesize, sd = sd_noise)
    data = cbind( x+e1 , y+e2 )

  } else if ( name == "normal") {

    theta = rnorm(samplesize, mean=0, sd=1)
    theta = sort(theta)
    x = rnorm(samplesize, sd=2)
    x = sort(x)
    y = rnorm(samplesize, sd=1)

    data = cbind(x, y)

  } else if ( name == "bananas") {

    phi = runif(samplesize, min=0, max=2*pi)
    phi = sort(phi)
    x1 = (10)*cos(phi)
    y1 = (10)*sin(phi)

    x2 = x1 - 2*(x1>0)
    y2 = y1 + 10*(x1>0)
    true_mani = cbind(x2, y2)

    e1 = rnorm(samplesize, sd = 2)
    e2 = rnorm(samplesize, sd = .5)
    data = cbind( x2+e1 , (y2+e2)/3 )

  } else if ( name == "line") {

    x = seq(from = -3, to = 3, length.out = samplesize)
    e1 = rnorm(40, sd=.0001)
    e2 = rnorm(40, sd=.0001)
    true_mani = cbind(x, rep(0,samplesize))

    data = cbind( x + e1, e2 )

  } else if ( name == "two-clouds") {

    x = rnorm(samplesize, sd=1)
    y = rnorm(samplesize, sd=1)
    z = sample(c(-3,3),samplesize,replace=T)

    data = cbind(x+z, y)

  } else if ( name == "not-dense-spiral") {

    t = seq(0,5*pi, length.out = samplesize)
    x = 1 * t * cos(t)
    y = 1 * t * sin(t)
    true_mani = cbind(x, y)

    e1 = rnorm(length(t), sd = sd_noise)
    e2 = rnorm(length(t), sd = sd_noise)
    data = cbind( x+e1 , y+e2 )

  } else if ( name == "sin-cos-curve") {

    t = seq(-pi/2,pi/2, length.out = samplesize)
    x = 2*sin(t)
    y = cos(t)+cos(2*t)+cos(6*t)
    true_mani = cbind(x, y)

    e1 = rnorm(length(t), sd = sd_noise)
    e2 = rnorm(length(t), sd = sd_noise)
    data = cbind( x+e1 , y+e2 )

  } else if ( name == "hetero-spiral") {

    t = seq(0,5*pi, by=.1)
    x = 1 * t^(3/2) * cos(t)
    y = 1 * t^(3/2) * sin(t)
    true_mani = cbind(x, y)

    e1 = rnorm(length(t), sd=t^(3/2)/10)
    e2 = rnorm(length(t), sd=t^(3/2)/10)
    data = cbind( x+e1 , y+e2 )

  } else if ( name == "two-peaks") {

    x = rnorm(samplesize, sd=1)
    y = rnorm(samplesize, sd=0.5)
    z = sample(c(-2,2),samplesize,replace=T)

    data = cbind(x+z, y)

  }
  if(plotTrue){
    plot(true_mani, pch=19, xlab='', ylab='', main = paste("true manifold", name, sep = " "))
    plot(data, pch=19, xlab='', ylab='', main = paste("noisy manifold", name, sep = " "))
  }

  result = list(data = data, true_mani = true_mani)
  return(result)

}













