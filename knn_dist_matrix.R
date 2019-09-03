




#############################################################################################
# knn provides with the classification with knn method for a given distance matrix
#############################################################################################

knn_dist_matrix <- function(ind_train,group,ksNN,lksNN, dist_matrix,N){
 
  ind_train<- sort(ind_train)
  groups<- group[ind_train]
  ind_seq<- 1:N
  ind_test<- ind_seq[-ind_train]
  n <- length(ind_train)
  n.new<- length(ind_test)
  good.class.l2 <- matrix(NA,nrow=lksNN,ncol=1)
  
  kkNN <- 1
  for (k in ksNN){
    
    groups.l2 <- matrix(NA,nrow=n,ncol=1)
    
    for (i in 1 : n){
      ind_temp<- ind_train[i]
      dist.l2 <- matrix(dist_matrix[ind_temp,ind_train],nrow=n,ncol=1)
      dist.l2[i] <- max(dist.l2)
      
      dist.l2.s <- sort(dist.l2,index.return = TRUE)
      g1 <- sum(groups[dist.l2.s$ix[1:k]]==1)
      g2 <- sum(groups[dist.l2.s$ix[1:k]]==2)
      if (g1 > g2){groups.l2[i] <- 1}else{groups.l2[i] <- 2}
      
    }
    good.class.l2[kkNN] <- sum(groups.l2==groups)
    kkNN <- kkNN + 1
    
  }
  
  k.knn.l2 <- min(which(good.class.l2==max(good.class.l2)))
  groups.knn.l2 <- matrix(NA,nrow=n.new,ncol=1)
  
  for (i in 1 : n.new){
    ind_temp<- ind_test[i]
    dist.l2 <- matrix(dist_matrix[ind_temp,ind_train],nrow=n,ncol=1)

    dist.l2.s <- sort(dist.l2,index.return = TRUE)
    g1 <- sum(groups[dist.l2.s$ix[1:k.knn.l2]]==1)
    g2 <- sum(groups[dist.l2.s$ix[1:k.knn.l2]]==2)
    if (g1 > g2){groups.knn.l2[i] <- 1}else{groups.knn.l2[i] <- 2}
  }
  
  return(length(which(groups.knn.l2!=group[ind_test]))/n.new)
  
}

knn_dist_matrix_k_optim<- function(ind_train,group,ksNN,lksNN, dist_matrix,N){
  
  ind_train<- sort(ind_train)
  groups<- group[ind_train]
  ind_seq<- 1:N
  ind_test<- ind_seq[-ind_train]
  n <- length(ind_train)
  n.new<- length(ind_test)
  good.class.l2 <- matrix(NA,nrow=lksNN,ncol=1)
  
  kkNN <- 1
  for (k in ksNN){
    
    groups.l2 <- matrix(NA,nrow=n,ncol=1)
    
    for (i in 1 : n){
      ind_temp<- ind_train[i]
      dist.l2 <- matrix(dist_matrix[ind_temp,ind_train],nrow=n,ncol=1)
      dist.l2[i] <- max(dist.l2)
      
      dist.l2.s <- sort(dist.l2,index.return = TRUE)
      g1 <- sum(groups[dist.l2.s$ix[1:k]]==1)
      g2 <- sum(groups[dist.l2.s$ix[1:k]]==2)
      if (g1 > g2){groups.l2[i] <- 1}else{groups.l2[i] <- 2}
      
    }
    good.class.l2[kkNN] <- sum(groups.l2==groups)
    kkNN <- kkNN + 1
    
  }
  
  k.knn.l2 <- min(which(good.class.l2==max(good.class.l2)))
  groups.knn.l2 <- matrix(NA,nrow=n.new,ncol=1)
  
  for (i in 1 : n.new){
    ind_temp<- ind_test[i]
    dist.l2 <- matrix(dist_matrix[ind_temp,ind_train],nrow=n,ncol=1)
    
    dist.l2.s <- sort(dist.l2,index.return = TRUE)
    g1 <- sum(groups[dist.l2.s$ix[1:k]]==1)
    g2 <- sum(groups[dist.l2.s$ix[1:k]]==2)
    if (g1 > g2){groups.knn.l2[i] <- 1}else{groups.knn.l2[i] <- 2}
  }
  
  return(k.knn.l2)
  
}

