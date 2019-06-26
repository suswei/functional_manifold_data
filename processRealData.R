
library(fda) # for growth dataset
library("fda.usc") # for tecator dataset

# returns a list of lists where every element describes a different dataset. all datasets are smoothed

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
  detach(growth)
  
  ### Tecator from fda package
  # This is from Mahalanobis technometrics paper
  # "The classification problem here is to separate meat samples with a high fat content(more than 20%) from samples with a low fat content (less than20%) based on absorbance. 
  # Among the 215 samples, 77 have high fat content and 138 have low fat content."
  
  
  # Warning: there are some weird boundary effects here!
  # Mahalanobis paper regarding smoothing
  # "second-order derivatives of the observed functions produces lower misclassification rates. 
  # Therefore, the analysis of the original data and their second-order derivatives are carried out. 
  # In both cases, the original discrete observations and their second differences are converted to 
  # functional observations using a B- spline basis of order 6 with 20 and 40 basis functions, respectively."
  
  data(tecator)
  absorp <- tecator$absorp.fdata
  Fat20 <- ifelse(tecator$y$Fat < 20, FALSE, TRUE) # class membership
  
  # bspline parameters
  norder <- 6
  lambda <- 1e-2
  
  
  for (nbasis in c(20,40)) {
    basis <- create.bspline.basis(tecator[["absorp.fdata"]]$rangeval, nbasis, norder)
    growfdPar <- fdPar(basis, 0, lambda)
    
    tecator_fd <- smooth.basis(tecator[["absorp.fdata"]]$argvals, t(absorp$data), growfdPar)$fd
    grid = absorp$argvals
    # grid_chopped = grid[-c(1,length(grid))] #TODO: is this cheating? we have to do it to avoid boundary effects
    if (nbasis==20){
      # Evaluate on original time grid
      list_of_lists[[list_index]] = list(data = eval.fd(grid, tecator_fd), reg_grid = grid, true_group = Fat20, name = "tecator")
      list_index = list_index+1
    } else{
      # Evaluate on original time grid
      list_of_lists[[list_index]] = list(data = eval.fd(grid, tecator_fd, 2), reg_grid = grid, true_group = Fat20, name = "tecator_der2")
      list_index = list_index+1
    }
  }
  
  # Lin and Yao's contamination paper regarding tecator data preporcessing:
  # we predict the fat content based on the first derivative curves approximated 
  # by the difference quotient between measurements at adjacent wavelengths, shown in the left panel of Figure 2.
  tecatorLY_fd=fdata.deriv(absorp,nderiv=1,method="diff",class.out='fdata')
  list_of_lists[[list_index]] = list(data = t(tecatorLY_fd$data), reg_grid = tecatorLY_fd$argvals, true_group = Fat20, name = "tecator_LYdiff1storder")
  list_index = list_index+1
  
  return(list_of_lists)

}