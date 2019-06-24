library("rlist")
processRealData = function(){
  
  list_of_lists = list()
  
  ### GROWTH CURVES from fda package
  # A list containing the heights of 39 boys and 54 girls from age 1 to 18 and the ages at which they were collected.
  
  #  --- Smooth, penalizing the Nth derivative (Lfdobj)  --
  #  A standard non-monotone smoothing using function smooth.basis
  #      with a penalty on the Nth derivative (Lfdobj) in order to estimate a
  #      smooth acceleration function
  #  This gives a smoother estimate of the acceleration functions
  #  Smoothing code follows https://github.com/cran/fda/blob/master/demo/growth.R with modification to simultaneous smoothing of both genders
  
  attach(growth)
  true_group<- c(rep(TRUE,39),rep(FALSE,54)) # TRUE = boy, FALSE = girl
  (ageRng <- range(age))
  agefine <- seq(ageRng[1],ageRng[2],length=101)
  knots  <- age
  norder <- 6
  nbasis <- length(knots) + norder - 2
  hgtbasis <- create.bspline.basis(range(knots), nbasis, norder, knots)
  
  list_index = 1
  for (Lfdobj in 1:4) {
    lambda <- 1e-2
    growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)
    
    # Need 'hgtm', 'hgtf', e.g., from attach(growth)
    hgtfd <- smooth.basis(growth$age, cbind(growth$hgtm,growth$hgtf), growfdPar)$fd
    growthhat <- eval.fd(agefine, hgtfd) # growth 
    velmhat <- eval.fd(agefine, hgtfd, 1) # velocity
    accmhat <- eval.fd(agefine, hgtfd, 2) # acceleration
    
    growthhat_list = list(data = growthhat, reg_grid = agefine, true_group = true_group, name = sprintf("growthhat_der%d",Lfdobj))
    list_of_lists[[list_index]] = growthhat_list
    velmhat_list = list(data = velmhat, reg_grid = agefine, true_group = true_group, name = sprintf("velmhat_der%d",Lfdobj))
    list_of_lists[[list_index+1]] = velmhat_list
    accmhat_list = list(data = accmhat, reg_grid = agefine, true_group = true_group, name = sprintf("accmhat_der%d",Lfdobj))
    list_of_lists[[list_index+2]]= accmhat_list
    list_index = list_index + 3
  }
  
  return(list_of_lists)
}