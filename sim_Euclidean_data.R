## Input 
# name of the Euclidean manifold
# samplesize
# SNR 
# reg_sampling : TRUE = regular sampling, FALSE = random uniform
# plotTrue  

## Output
# true_data : samplesize x K matrix containing the original data (no noise)
# discrete_data : samplesize x K matrix containing the observed data
# true_geo : samplesize x samplesize matrix containing true pairwise geo disctance

#Example of call
# data<- sim_Euclidean_data('circle',100,25,TRUE,TRUE)

sim_Euclidean_data<-function(name,samplesize,SNR,reg_sampling,plotTrue){
  obj = EuclideanExamples(name, samplesize, SNR, plotTrue,reg_sampling,300)
  data = data.matrix(obj$data)
  true_mani  = data.matrix(obj$true_mani)
  true_geo = data.matrix(obj$true_geo)
  return(list('true_data'=true_mani,'discrete_data'=data,'true_geo'=true_geo))
}


