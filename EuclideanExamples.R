# generate noisy manifold data
# TODO: not sure if all of these noises are "normal" to the manifold
# TODO: need to add the true manifold to all the examples, have only done this for circle
# TODO: clean up extra parameters that don't get called when plotting is turned off
# options:
#   spiral
#   cross
#   circle
#   right-angle
#   half-moon
#   normal
#   bananas

EuclideanExamples <- function(name,samplesize){

  if ( name == "spiral") {

    h= 1.4
    n.grid = 2

    n.steps = samplesize
    t = seq(pi,5*pi, length.out = n.steps)
    x = 1 * t * cos(t)
    y = 1 * t * sin(t)
    true_mani = data.frame(x = x, y = y)

    # TODO: is this noise normal to the manifold?
    e1 = rnorm(length(t), sd=.3)
    e2 = rnorm(length(t), sd=.3)
    data = matrix(NA, nrow=n.steps, ncol=n.grid)
    data[,1] = x+e1
    data[,2] = y+e2

    length.lm = 120/h^2
    d=1
    grid = 1
    x.from = -21
    x.to = 18
    y.from = -15
    y.to = 25

  } else if ( name == "manifold") {

    # calling manifold_data.py will return double the samplesize,
    # first set is the noisy manifold data, second set is the true manifold
    source_python('manifold_data.py')
    stacked = manifold_data(samplesize)

    data = stacked[1:samplesize, ]
    true_mani = stacked[-(1:samplesize), ]

    x.from = -2
    x.to = 16
    y.from = -13
    y.to = 3
    grid = 1
    length.lm = 1


  } else if ( name == "cross") {

    x = rnorm(samplesize, sd=2)
    y = rnorm(samplesize, sd=.2)

    data = data.frame(x=x,y=y)

    x = rnorm(samplesize, sd=.2)
    y = rnorm(samplesize, sd=2)

    data = rbind(data, data.frame(x=x,y=y))

    h= 1
    length = 8/(h^2)
    length.lm = 8/(h^2)
    d=1
    grid = .5

    x.from = -20
    x.to = 20
    y.from = -10
    y.to = 10

  } else if ( name == "circle") {

    ###########################
    ##circle

    angle = runif(samplesize, 0,2*pi)
    true_mani = data.frame(x = 10*cos(angle), y = 10*sin(angle))

    # add noise normal to manifold
    alpha = rnorm(samplesize, sd=2)
    data = data.frame(x = (10+alpha)*cos(angle) , y= (10+alpha)*sin(angle) )

    h = 2
    length.lm = 200/h^2
    d=1
    h= 1
    grid = .8
    x.from = -15
    x.to = 15
    y.from = -15
    y.to = 15

  } else if ( name == "right-angle") {


    ###########################
    ##right-angle

    h= 2

    n.grid = 2
    n.steps = samplesize
    t = seq(0,5*pi, length.out = n.steps)

    x = c( seq(-10,0,length.out = n.steps/2), rep(0, n.steps/2)  )
    y = c( rep(0, n.steps/2), seq(0,10,length.out = n.steps/2) )


    e1 = rnorm(n.steps, sd=.1)
    e2 = rnorm(n.steps, sd=.1)

    data = matrix(NA, nrow=n.steps, ncol=n.grid)
    data[,1] = x+e1
    data[,2] = y+e2


    length.lm = 8/h^2
    d=1
    grid = .35

    x.from = -12
    x.to = 2
    y.from = -2
    y.to = 12

  } else if ( name == "half-moon") {

    ##############################
    ##half moon

    h=.25

    n.grid = 2
    n.steps = samplesize

    length.lm = 1/h^2
    d=1
    grid = .1

    x.from = -1.2
    x.to = 1.2
    y.from = -1.2
    y.to = 1.2

    theta = rnorm(n.steps, mean=0, sd=1)
    theta = sort(theta)
    x = cos(theta)
    y = sin(theta)

    e1 = rnorm(n.steps, sd=.1)
    e2 = rnorm(n.steps, sd=.1)
    data = matrix(NA, nrow=n.steps, ncol=n.grid)
    data[,1] = x+e1
    data[,2] = y+e2

    center = c(1, 0)
    sigma.append = .07

  } else if ( name == "normal") {

    #################################33
    ##normal distribution



    h=1.5

    n.grid = 2
    n.steps = samplesize

    length.lm = 3/h^2
    d=1
    grid = .5

    x.from = -8
    x.to = 8
    y.from = -4
    y.to = 4

    theta = rnorm(n.steps, mean=0, sd=1)
    theta = sort(theta)
    x = rnorm(n.steps, sd=2)
    x = sort(x)
    y = rnorm(n.steps, sd=1)

    data = matrix(c(x,y), nrow=n.steps, ncol=n.grid)

    center = c(0, 0)
    sigma.append = .25

  } else if ( name == "bananas") {

    #############################
    ##two bananas

    h= 1

    n.steps = samplesize
    phi = runif(n.steps, min=0, max=2*pi)
    phi = sort(phi)

    x1 = (10)*cos(phi)
    y1 = (10)*sin(phi)

    x2 = x1 - 2*(x1>0)
    y2 = y1 + 10*(x1>0)

    e1 = rnorm(n.steps, sd=2)
    e2 = rnorm(n.steps, sd=.5)

    data = matrix( c( x2+e1 , y= (y2+e2)/3), ncol=2 )

    length.lm = 15/h^2

    d=1
    grid = 1
    x.from = -18
    x.to = 15
    y.from = -15
    y.to = 15

    center = c(-10, 0)
    sigma.append = .5

  } else if ( name == "line") {

    ###############
    ##line

    x = seq(from = -3, to = 3, length.out = samplesize)
    e1 = rnorm(40, sd=.0001)
    e2 = rnorm(40, sd=.0001)

    data = data.frame(x = x +e1, y=0+e2 )

    length = .5
    bf.length = .3
    length.lm = bf.length
    d=1
    h= .8
    grid = .1

    x.from = -4
    x.to = 4
    y.from = -3
    y.to = 3

  } else if ( name == "two-clouds") {
    ####################
    ##two clouds
    x = rnorm(samplesize, sd=1)
    y = rnorm(samplesize, sd=1)
    z = sample(c(-3,3),samplesize,replace=T)

    data = data.frame(x=x+z,y=y)

    length = 2
    length.lm = length;

    d=1
    h= 2
    grid = 0.4

    x.from = -6
    x.to = 6
    y.from = -4
    y.to = 4


  } else if ( name == "Archimedean-spiral") {
    ########################
    # Archimedean spiral


    t = seq(0,5*pi, length.out = samplesize)

    x = 1 * t * cos(t)
    y = 1 * t * sin(t)


    e1 = rnorm(length(t), sd=.3)
    e2 = rnorm(length(t), sd=.3)

    data = data.frame(x = x+e1 , y= y+e2 )

    length = 10
    bf.length = .2
    length.lm = bf.length
    d=1
    h= 1.5
    grid = .5

    x.from = -21
    x.to = 18
    y.from = -15
    y.to = 25

  } else if ( name == "not-dense-spiral") {
    ##############################################
    # spiral not dense

    t = seq(0,5*pi, length.out = samplesize)

    x = 1 * t * cos(t)
    y = 1 * t * sin(t)


    e1 = rnorm(length(t), sd=.3)
    e2 = rnorm(length(t), sd=.3)

    data = data.frame(x = x+e1 , y= y+e2 )

    length = 10
    bf.length = .2
    length.lm = bf.length
    d=1
    h= 1.2
    grid = .5

    x.from = -21
    x.to = 18
    y.from = -15
    y.to = 25

  } else if ( name == "sin-cos-curve") {

    #############################
    ##  sin-cos curve

    t = seq(-pi/2,pi/2, length.out = samplesize)

    x = 2*sin(t)
    y = cos(t)+cos(2*t)+cos(6*t)

    e1 = rnorm(length(t), sd=.01)
    e2 = rnorm(length(t), sd=.01)

    data = data.frame(x = x+e1 , y= y+e2 )
    plot(y~x, data)

    length = 14
    unit.length = .05
    length.lm = unit.length
    d=1
    h= .05
    grid = .11

    x.from = -3
    x.to = 3
    y.from = -2.25
    y.to = 3.25

  } else if ( name == "hetero-spiral") {

    ###########################
    # heteroscedastic spiral

    t = seq(0,5*pi, by=.1)
    x = 1 * t^(3/2) * cos(t)
    y = 1 * t^(3/2) * sin(t)


    e1 = rnorm(length(t), sd=t^(3/2)/10)
    e2 = rnorm(length(t), sd=t^(3/2)/10)

    data = data.frame(x = x+e1 , y= y+e2 )

    length = 10
    unit.length = 1
    length.lm = unit.length
    d=1
    h= 20
    grid = 2

    x.from = -80
    x.to = 60
    y.from = -50
    y.to = 70

  } else if ( name == "two-peaks") {

    ########################################
    # two peaks

    x = rnorm(samplesize, sd=1)
    y = rnorm(samplesize, sd=0.5)
    z = sample(c(-2,2),samplesize,replace=T)

    data = data.frame(x=x+z,y=y)

    length = 2
    unit.length = 0.1
    length.lm = unit.length
    d=1
    h= 0.3
    grid = 0.2

    x.from = -6
    x.to = 6
    y.from = -15
    y.to = 15
  }

  plot(true_mani, pch=19, xlab='', ylab='', main = paste("true manifold", name, sep = " "))
  plot(data, pch=19, xlab='', ylab='', main = paste("noisy manifold", name, sep = " "))

  # length.lm only used for plotEigenvectorField
  result = list(data = data, true_mani = true_mani,
                x.from = x.from, x.to = x.to, y.from = y.from,
                y.to = y.to, grid = grid,
                length.lm = length.lm, name = name)
  return(result)
}













