##### Function that assess how good a particular estimation of a geodesic matrix is good

# for the moment we just use the relative error || hat A - A ||_F / || A ||_F

assess_goodness_estimation <- function(estim_mat,true_geo){

  rmse = sqrt(sum((estim_mat-true_geo)^2)/sum(true_geo^2))

  epsilons = seq(1,100,1)
  epsilon_isometry = lapply(epsilons,check_epsilon_isometry,estim_mat=estim_mat,true_geo=true_geo)

  smallest_epsilon = min(epsilons[unlist(epsilon_isometry)]) #find out which is the smallest epsilon such that isometry holds, returns Inf if nothing in the eipsilons exhibit this behavior

  return(list("rmse" = rmse, "smallest_epsilon" = smallest_epsilon))
}

check_epsilon_isometry <- function(estim_mat,true_geo,epsilon){

  all( ( (1-epsilon)*true_geo <= estim_mat) & (estim_mat <= (1+epsilon) *true_geo) )

}
