##### Function that assess how good a particular estimation of a geodesic matrix is good

# for the moment we just use the relative error || hat A - A ||_F / || A ||_F

assess_goodness_estimation <- function(estim_mat,true_geo){
  sqrt(sum((estim_mat-true_geo)^2)/sum(true_geo^2))
}

