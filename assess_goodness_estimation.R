##### Function that assess how good a particular estimation of a geodesic matrix is good

# for the moment we just use the relative error || hat A - A ||_F / || A ||_F

assess_goodness_estimation <- function(estim_mat,true_geo){

  rmse = sqrt(sum((estim_mat-true_geo)^2)/sum(true_geo^2))

  epsilons = seq(0,1,0.01)
  epsilon_isometry_prop = lapply(epsilons,check_epsilon_isometry,estim_mat=estim_mat,true_geo=true_geo)
  epsilon_isometry_auc = auc(epsilons, unlist(epsilon_isometry_prop), type = "spline")
  
pearson_corr = cor.test(estim_mat[lower.tri(estim_mat)],true_geo[lower.tri(true_geo)], method="pearson")$estimate

  return(list("rmse" = rmse, "epsilon_isometry_auc" = epsilon_isometry_auc, "pearson_corr"=pearson_corr))
}

check_epsilon_isometry <- function(estim_mat,true_geo,epsilon){

  mean( ( (1-epsilon)*true_geo <= estim_mat) & (estim_mat <= (1+epsilon) *true_geo) )

}
