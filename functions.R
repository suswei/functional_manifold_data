################################################
## calculates the local covariance at (x0,y0)
loc.cov = function(obj,x0,h,shift_type,kernel){

  data = data.matrix(obj$data)

  weights = exp( -((data[,1]-x0[1])^2+(data[,2]-x0[2])^2)/h^2 ) / nrow(data)

  if (shift_type == "ms") {

    local.mean.x = sum(weights * data[,1]) / sum(weights)
    local.mean.y = sum(weights * data[,2]) / sum(weights)

  } else if(shift_type == "scms") {

    kernel_sigma = h
    kernelscms = scms$make_isotropic_gaussian_kernel(kernel_sigma)
    denoised = scms$subspace_constrained_mean_shift_update(x0, data, kernelscms, kernel_sigma)
    local.mean.x = x0[1]+denoised[1]
    local.mean.y = x0[2]+denoised[2]

  } else if(shift_type == "none") {
    local.mean.x = x0[1]; local.mean.y = x0[2];
  }

  if (kernel == "double-peak") {
    # double-peak kernel
    m = matrix(c(sum(weights*(data[,1]-local.mean.x)^2) ,
                 sum(weights*(data[,1]-local.mean.x)*(data[,2]-local.mean.y)) ,
                 sum(weights*(data[,2]-local.mean.y)*(data[,1]-local.mean.x)) ,
                 sum(weights*(data[,2]-local.mean.y)^2  ) ) ,nrow=2)
  } else if (kernel == "butterfly") {
    # butterfly kernel
    m = matrix(c(sum(weights*(data[,1]-local.mean.x)^2 / ( (data[,1]-local.mean.x)^2+(data[,2]-local.mean.y)^2 ) ) ,
                 sum(weights*(data[,1]-local.mean.x)*(data[,2]-local.mean.y) / ( (data[,1]-local.mean.x)^2+(data[,2]-local.mean.y)^2 ) ) ,
                 sum(weights*(data[,2]-local.mean.y)*(data[,1]-local.mean.x) / ( (data[,1]-local.mean.x)^2+(data[,2]-local.mean.y)^2 ) ) ,
                 sum(weights*(data[,2]-local.mean.y)^2 / ( (data[,1]-local.mean.x)^2+(data[,2]-local.mean.y)^2 ) )) ,nrow=2)
  }
  return(m)
}


# plot eigenvector field
plotEigenvectorField = function(obj,h,shift_type,kernel) {

  plot(obj$data, pch=19, xlab='', ylab='',main=obj$name)

  for (x0 in seq(from=obj$x.from,to=obj$x.to,by=obj$grid)){
    for (y0 in seq(from=obj$y.from,to=obj$y.to,by=obj$grid)){

      ei = eigen(loc.cov(obj,c(x0,y0),h,shift_type,kernel))

      segments(x0-ei$vectors[1,2]*ei$values[2]*obj$length.lm,
               y0-ei$vectors[2,2]*ei$values[2]*obj$length.lm,
               x0+ei$vectors[1,2]*ei$values[2]*obj$length.lm,
               y0+ei$vectors[2,2]*ei$values[2]*obj$length.lm,
               col="red", lwd=3)
      segments(x0-ei$vectors[1,1]*ei$values[1]*obj$length.lm,
               y0-ei$vectors[2,1]*ei$values[1]*obj$length.lm,
               x0+ei$vectors[1,1]*ei$values[1]*obj$length.lm,
               y0+ei$vectors[2,1]*ei$values[1]*obj$length.lm,
               col="blue", lwd=3)
    }
  }
}

anglebetween = function(a,b) {
  anglebetween = acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
}

eigenvectorFlow = function(obj,x0,shift_type,kernel,h,accelerated) {

  #### follow

  lambda = 1
  n = 100000

  path = array(NA, dim=c(2,n,2))
  path[1,1,] = x0
  path[2,1,] = x0

  # j=1 forward, j=2 backward
  for (j in 1:2){
    # we need to do the first step explicitely to set a direction of travel
    tensor = loc.cov(obj,path[j,1,],h,shift_type,kernel)
    eig = eigen(tensor)
    direction = eig$vectors[,1]

    step.size = eig$values[1] * lambda
    prev.step = direction * step.size

    if (j==2) {
      prev.step = -prev.step
    }

    path[j,2,] = path[j,1,] + prev.step #this step is the same for vanilla and accelerated

    # cycle for all other steps
    #for (i in 3:n){
    i = 3
    while (sqrt(sum(prev.step^2)) > 0.001 & i<(n-1)) {


      if (accelerated){
        y  = path[j,i-1,] + (i-2)/(i+1) * (path[j,i-1,] - path[j,i-2,])

        htemp = 2
        data = data.matrix(obj$data)
        weights = exp( -((data[,1]-y[1])^2+(data[,2]-y[2])^2)/htemp^2 ) / nrow(data)
        y[1] = sum(weights * data[,1]) / sum(weights)
        y[2] = sum(weights * data[,2]) / sum(weights)

        tensor = loc.cov(obj,y,h,shift_type,kernel)
        eig = eigen(tensor)
        new.step = eig$vectors[,1]
        step.size = eig$values[1] * lambda

        # if the angle obtuse, we change the direction
        if (new.step %*% prev.step < 0) { new.step = -new.step }

        step.size = 1/i
        path[j,i,] = y + step.size*new.step
        prev.step = new.step
      }
      else
      {
        tensor = loc.cov(obj,path[j,i-1,],h,shift_type,kernel)
        eig = eigen(tensor)
        direction = eig$vectors[,1]
        step.size = eig$values[1] * lambda
        new.step = direction * step.size
        # if the angle obtuse, we change the direction
        if (new.step %*% prev.step < 0) { new.step = -new.step }
        path[j,i,] = path[j,i-1,] + new.step
        prev.step = new.step
      }
      i = i+1
    }


  }

  data = obj$data

  # redraw plot
  plot(data, pch=19, xlab='', ylab='', main=paste("eigflow",'h=',h, sep=' ') )

  # draw point where I start
  points(x0[1],x0[2], col="green", pch=19, cex=1.5)

  # draw integral curves
  for (j in 1:2) segments( path[j, 1:(n-1),1], path[j, 1:(n-1),2], path[j,2:n,1], path[j,2:n,2], col=rainbow(n-1), lwd=3 )

  return(path)

}
