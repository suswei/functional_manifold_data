approx.fourier <- function(DATA, nbasis, Range.grid, period=NULL)
{
  #####################################################################
  # Fourier approximation of the curves containing in DATA :
  # -----------------------------------------------------------
  # "nbasis" and "period" allow to define the Fourier basis
  # "COEF[, i]" corresponds to the Fourier expansion
  # of the discretized curve contained in DATA[i, ]. 
  # The Fourier approximation of the curve contained in "DATA[i, ]" 
  # is given by "FOURIER %*% COEF[, i]"
  #####################################################################
  p <- ncol(DATA)
  nbasismax <- (p - nbasis)%/%2
  if(nbasis > nbasismax){
    stop(paste("give a number nbasis smaller than ",nbasismax, " for avoiding ill-conditioned matrix"))
  }
  a <- Range.grid[1]
  b <- Range.grid[2]
  x <- seq(a, b, length = p)
  if(is.null(period)) period <- b - a
  FOURIER <- fourier(x, nbasis, period)
  CMAT <- crossprod(FOURIER)
  DMAT <- crossprod(FOURIER, t(DATA))
  COEF <- symsolve(CMAT, DMAT)
  list(COEF = COEF, APPROX = FOURIER %*% COEF)
}
#####################################################################
#####################################################################
#####################################################################
approx.spline.deriv <- function(DATA, nderivs, nknot, Range.grid)
{
  #####################################################################
  # B-spline approximation of the successive derivatives of the curves 
  # containing in DATA :
  # -----------------------------------------------------------
  # "nderivs" defines the number of the wished derivative
  # "nknot" allows to define the B-spline basis
  # "COEF[, i]" corresponds to the B-spline expansion
  # of the discretized curve contained in DATA[i, ]. 
  # The B-spline approximation of the successive derivatives the curve 
  # contained in "DATA[i, ]" is given by "APPROX[, i]"
  #####################################################################
  library(splines)
  p <- ncol(DATA)
  a <- Range.grid[1]
  b <- Range.grid[2]
  Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
  x <- seq(a, b, length = p)
  ord <- nderivs + 3
  nknotmax <- (p - ord - 1)%/%2
  if(nknot > nknotmax){
    stop(paste("give a number nknot smaller than ",nknotmax, " for avoiding ill-conditioned matrix"))
  }
  Delta <- sort(c(rep(c(a, b), ord), Knot))
  BSPLINE <- splineDesign(Delta, x, ord)
  CMAT <- crossprod(BSPLINE)
  DMAT <- crossprod(BSPLINE, t(DATA))
  COEF <- symsolve(CMAT, DMAT)
  BSPLINEDER <- splineDesign(Delta, x, ord, rep(nderivs, length(x)))
  list(COEF = COEF, APPROX = BSPLINEDER %*% COEF)
}
#####################################################################
#####################################################################
#####################################################################
approx.spline <- function(DATA, order, nknot, Range.grid)
{
  #####################################################################
  # B-spline approximation of the curves containing in DATA :
  # -----------------------------------------------------------
  # "order" and "nknot" allow to define the B-spline basis
  # "COEF[, i]" corresponds to the B-spline expansion
  # of the discretized curve contained in DATA[i, ]. 
  # The B-spline approximation of the curve contained in "DATA[i, ]" 
  # is given by "APPROX[, i]"
  #####################################################################
  library(splines)
  nknotmax <- (p - order - 1)%/%2
  if(nknot > nknotmax){
    stop(paste("give a number nknot smaller than ",nknotmax, " for avoiding ill-conditioned matrix"))
  }
  p <- ncol(DATA)
  a <- Range.grid[1]
  b <- Range.grid[2]
  Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
  x <- seq(a, b, length = p)
  Delta <- sort(c(rep(c(a, b), order), Knot))
  BSPLINE <- splineDesign(Delta, x, order)
  CMAT <- crossprod(BSPLINE)
  DMAT <- crossprod(BSPLINE, t(DATA))
  COEF <- symsolve(CMAT, DMAT)
  list(COEF = COEF, APPROX = BSPLINE %*% COEF)
}
#####################################################################
#####################################################################
#####################################################################
classif.bw <- function(PROBCURVES, Bw.seq, Group)
{
  ###############################################################
  # Computes a bandwidth among a fixed sequence Bw.seq 
  # from probability curves contained in the matrix 
  # "PROBCURVES" which minimizes the entropy:
  #   "PROBCURVES" contains the probability curves (row by row)
  #   "Bw.seq"     contains a fixed sequence of bandwidths
  #   "Group"      gives the unit numbers of the concerned group
  # Returns a list containing:
  #   "$index" is the index correponding to the optimal bandwidth
  #   "$bw"    is the optimal bandwidth 
  ###############################################################
  number.of.bw <- length(Bw.seq)
  entro <- 0
  options(warn=-1)
  for(j in 1:number.of.bw){
    prob.points <- PROBCURVES[Group,j]
    pp.start <- min(prob.points)
    pp.end <- max(prob.points)
    if(pp.start<pp.end){
      est <- density(prob.points, bw="bcv",from=pp.start, to=pp.end)
      entro[j] <- entropy(est$y, est$x)
    } else {
      entro[j] <- 0
    }
  }
  index <- order(entro[entro>0])[1]+sum(entro==0)
  return(list(index=index, bw=Bw.seq[index]))
}
#####################################################################
#####################################################################
#####################################################################
classif.hi <- function(CURVES, SEMIMETRIC, kind.of.kernel, bw, Group, semimetric, centrality,...)
{
  #####################################################################
  # Returns the "Heterogeneity Index (HI)" of the class
  # defined by the vector "Group" of the curves dataset "CURVES":
  #   "CURVES" matrix containing the curves dataset (row by row)
  #   "SEMIMETRIC" matrix containing the proximities between curves
  #   "bw"     is a selected bandwidth (see function "classif.bw")
  #   "Group"  gives the unit numbers of the concerned class
  #   "semimetric" gives the semimetric type required
  #            (see functions "semimetric.deriv", "semimetric.pca",...)
  #   "centrality": character string giving the centrality feature 
  #            compared with the mode (centrality="mean" or "median").
  #   "..."    arguments for defining the used semimetric
  #            computing HI
  # Returns the value of HI for the class defined by "Group"
  #####################################################################
  lg <- length(Group)
  p <- ncol(CURVES)
  rank.mode <- funopa.mode(bw, SEMIMETRIC[Group, Group],kind.of.kernel)
  Modal.curve <- CURVES[Group[rank.mode],]
  if(centrality=="mean")
    Centrality.feature <- apply(CURVES[Group, ], 2, sum)/lg
  if(centrality=="median")
    Centrality.feature <- median.npfda(SEMIMETRIC[Group, Group])
  MM0 <- rbind(Modal.curve, Centrality.feature, rep(0, p))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  SM.MM0 <- sm(MM0, MM0, ...)
  return(SM.MM0[1,2]/(SM.MM0[1,3]+SM.MM0[2,3]))
}
#####################################################################
#####################################################################
#####################################################################
classif.npfda <- function(CURVES, ..., kind.of.kernel="quadratic", semimetric="deriv", threshold=0.1, nb.bw=100, nss=0, mspg=10, centrality="mean")
{
  #################################################################
  # Performs the unsupervised classification procedure:
  # "CURVES" contains the curves (row by row)
  # ... parameter needed for the used semimetric routine
  # "kind.of.kernel" gives the kernel used for computing the modes
  # "semimetric" contains the semimetric type (default="deriv")
  # "threshold"  contains the value which gives the minimum gain
  #    accepted in terms of splitting score
  # "nb.bw"      gives the number of bandwidths for building the
  #    corresponding sequence (default="100")
  # "nss"        gives the number of subsamples required for
  #    computing the subsampling heterogeneity index (default="50")
  #    (if "nss=0", the method uses HI instead of SHI which reduces
  #    the computational cost)
  # "mspg"       gives the minimal size per group authorized 
  # "centrality": character string giving the centrality feature 
  #     compared with the mode (centrality="mean" or "median") 
  #     in order to compute the heterogeneity index (default="mean")      
  # Returns a list:
  #  MEANS     matrix containing the mean curve of each group
  #  MEDIANS   matrix containing the median curve of each group
  #  MODES     matrix containing the modal curve of each group
  #  Bw.opt    vector containing the optimal bandwidth of each group
  #  Partition list giving the unit numbers componed each group
  #  Labels    character vector giving the classification tree
  #  Ssc       vector giving the splitting score of each accepted split 
  #####################################################################
  n <- nrow(CURVES)
  ##
  # Computes the proximity between curves according to the selected semimetric
  ############################################################################
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  SEMIMETRIC <- sm(CURVES, CURVES, ...)
  ##
  # Build the sequence of bandwidths
  ##################################
  Hrange <- range(SEMIMETRIC)
  Bw.seq <- seq(Hrange[1], Hrange[2] * 0.5, length = nb.bw + 1)[-1]
  ##
  # Probability curves of spectrometric data
  ##########################################
  PROBCURVES <- matrix(0, n, nb.bw)
  for(i in 1:n){
    PROBCURVES[i,] <- prob.curve(i, SEMIMETRIC, Bw.seq)
  }
  ##
  # Selecting the optimal bandwidth of full sample
  ################################################
  res.classif.bw <- classif.bw(PROBCURVES, Bw.seq, 1:n)
  Bw.opt <- res.classif.bw$bw
  index <- res.classif.bw$index
  ##
  # Computing SHI/HI of full sample
  #################################
  if(nss==0){
    shi <- classif.hi(CURVES, SEMIMETRIC, kind.of.kernel, Bw.opt, 1:n, semimetric, centrality, ...)
  } else {
    shi <- classif.shi(CURVES, SEMIMETRIC, kind.of.kernel, Bw.opt, 1:n, nss, semimetric, centrality, ...)
  }
  ##
  # Partitioning of full sample
  #############################
  current.part <- classif.part(CURVES, PROBCURVES, SEMIMETRIC, kind.of.kernel, index, Bw.seq, 1:n, shi, semimetric, threshold, nss, mspg, centrality, ...)
  current.partition <- list()
  current.partition[[1]] <- current.part
  current.part <- current.partition
  nbpart <- length(current.part)
  final.partition <- list()
  nbgroups <- length(current.part[[nbpart]]$Groups)
  ##
  # Recursive partitioning procedure
  ##################################
  if(nbgroups>0){
    Gains <- numeric(0)
    count1 <- 0
    Names <- list("")
    Names.new <- list("")
    Labels <- character(0)
    while(nbpart!=0){
      count3 <- 0
      test <- F
      for(l in 1:nbpart){
        count2 <- 0
        new.partition <- list()
        nbgroups <- length(current.part[[l]]$Groups)
        if(nbgroups>0){
          Logic3 <- F
          for(k in 1:nbgroups){
            label <- paste(Names[[l]],k,sep="")
            index <- current.part[[l]]$Index[k]
            Group <- current.part[[l]]$Groups[[k]]
            shi <- current.part[[l]]$SHI[k]
            new.partition[[k]] <- classif.part(CURVES, PROBCURVES, SEMIMETRIC, kind.of.kernel, index, Bw.seq, Group, shi, semimetric, threshold, nss, mspg, centrality, ...)
            logic1 <- (length(new.partition[[k]]$Groups)==0)
            if(!logic1){ logic2 <- (new.partition[[k]]$Ssc<threshold)
            } else {logic2 <- FALSE}
            if(logic1 || logic2){
              count1 <- count1 + 1
              Labels[count1] <- label
              final.partition[[count1]] <- current.part[[l]]$Groups[[k]]
              Bw.opt[count1] <- current.part[[l]]$Bw.opt[k]
            } else {
              test <- T
              count2 <- count2 + 1
              Names.new[[l]][count2] <- label
              count3 <- count3 + 1
              current.partition[[count3]] <- new.partition[[k]]
            }
          }
          Gains <- c(Gains, current.part[[l]]$Ssc)
          if(test){Names <- as.list(Names.new[[l]])}
        } else {
          count1 <- count1 + 1
          final.partition[[count1]] <- current.part[[l]]$Groups
        }
      }
      if(test){
        current.part <- current.partition
        nbpart <- length(current.part)
      } else {
        nbpart <- 0
      }
    }
    nb.modes <- length(final.partition)
    MODES <- matrix(0,nb.modes,ncol(CURVES))
    MEDIANS <- matrix(0,nb.modes,ncol(CURVES))
    MEANS <- matrix(0,nb.modes,ncol(CURVES))
    for(j in 1:nb.modes){
      Group <- final.partition[[j]]
      rank.mode <- funopa.mode(Bw.opt[j], SEMIMETRIC[Group, Group], kind.of.kernel)
      rank.median <- median.npfda(SEMIMETRIC[Group, Group])
      MODES[j,] <- CURVES[Group[rank.mode],]
      MEDIANS[j,] <- CURVES[Group[rank.median],]
      MEANS[j,] <- apply(CURVES[Group,], 2, sum)/length(Group)
    }
    return(list(MEANS=MEANS, MEDIANS=MEDIANS, MODES=MODES, Bw.opt=Bw.opt, Partition=final.partition, Labels=Labels, Ssc=Gains))
  } else { 
    cat("Classification procedure don't split the initial sample according to the selected criterion")
  }
}
#####################################################################
#####################################################################
#####################################################################
classif.part <- function(CURVES, PROBCURVES, SEMIMETRIC, kind.of.kernel, index, Bw.seq, Group, shi, semimetric, threshold, nss, mspg, centrality, ...)
{
  #################################################################
  # From the set of curves contained in the matrix "CURVES", it performs
  # a partition and computes several features
  #   "CURVES" contains the curves (row by row)
  #   "PROBCURVES" contains the probability curves (row by row)
  #   "SEMIMETRIC" matrix containing the proximities between curves
  #   ... parameter needed for the used semimetric routine
  #   "kind.of.kernel" gives the kernel used for computing the modes
  #   "index" is the index correponding to the optimal bandwidth
  #           (see function "classif.bw")
  #   "Bw.seq" contains a fixed sequence of bandwidths
  #   "Group"  gives the unit numbers of the concerned class
  #   "shi" subsampling heterogeneity index of the set of curves
  #         defined by the vector "Group"
  #   "semimetric" contains the semimetric type (default="deriv")
  #   "threshold"  contains the value which gives the minimum gain
  #      accepted in terms of splitting score
  #   "nss"        gives the number of subsamples required for
  #      computing the subsampling heterogeneity index (default="50")
  #      (if "nss=0", the method uses HI instead of SHI which reduces
  #      the computational cost)
  #   "mspg"       gives the minimal size per group authorized 
  #   "centrality": character string giving the centrality feature 
  #       compared with the mode (centrality="mean" or "median") 
  #       in order to compute the heterogeneity index (default="mean")      
  # Returns a list:
  #   "$Groups" list giving the unit numbers componed each group of the 
  #             performed partition
  #   "Bw.opt"  vector containing the optimal bandwidth of each group
  #   "$Index"  vector containing the index correponding to the optimal 
  #             bandwidth
  #   "$SHI"    vector giving the SHI for each group
  #   "$Ssc"    scalar giving the splitting score of performed partition
  #####################################################################
  res.classif.split <- classif.split(PROBCURVES, index, Group)
  nbgroups <- length(res.classif.split)
  res <- list(Groups=list())
  if(nbgroups!=0){
    minsize <- min(sapply(res.classif.split, length))
    if(minsize>mspg){
      Bw.opt <- 0
      Index <- 0
      shisub <- 0
      for(k in 1:nbgroups){
        Subgroup <- res.classif.split[[k]]
        res.classif.bw <- classif.bw(PROBCURVES, Bw.seq, Subgroup)
        Bw.opt[k] <- res.classif.bw$bw
        Index[k] <- res.classif.bw$index
        if(nss==0){
          shisub[k] <- classif.hi(CURVES, SEMIMETRIC, kind.of.kernel, Bw.opt[k], Subgroup, semimetric, centrality, ...)
        } else {
          shisub[k] <- classif.shi(CURVES, SEMIMETRIC, kind.of.kernel, Bw.opt[k], Subgroup, nss, semimetric, centrality, ...)  
        }
      }
      split.score <- 1- sum(sapply(res.classif.split, length)*shisub)/(length(Group)*shi)
      if(split.score>threshold){
        res <- list(Groups=res.classif.split, Bw.opt=Bw.opt, Index=Index, SHI=shisub, Ssc=split.score)
      }
    }
  }
  return(res)
}
#####################################################################
#####################################################################
#####################################################################
classif.shi <- function(CURVES, SEMIMETRIC, kind.of.kernel, bw, Group, nss, semimetric, centrality,...)
{
  #####################################################################
  # Returns the "Subsampling Heterogeneity Index (SHI)" of the class
  # defined by the vector "Group" of the curves dataset "CURVES":
  #   "CURVES" matrix containing the curves dataset (row by row)
  #   "SEMIMETRIC" matrix containing the proximities between curves
  #   "bw"     is a selected bandwidth (see function "classif.bw")
  #   "Group"  gives the unit numbers of the concerned class
  #   "nss"    is the fixed Number of SubSamples created for 
  #            computing SHI
  #   "semimetric" gives the semimetric type required
  #            (see functions "semimetric.deriv", "semimetric.pca",...)
  #   "centrality": character string giving the centrality feature 
  #            compared with the mode (centrality="mean" or "median").
  #   "..."    arguments for defining the used semimetric
  # Returns the value of SHI for the class defined by "Group"
  #####################################################################
  lg <- length(Group)
  size <- ceiling(lg/2)
  Hetero.index <- 0
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  p <- ncol(CURVES)
  for(e in 1:nss){
    Sample <- sort(sample(Group,size))
    SM.SAMPLE <- SEMIMETRIC[Sample, Sample]
    rank.mode <- funopa.mode(bw, SM.SAMPLE, kind.of.kernel)
    Modal.curve <- CURVES[Sample[rank.mode],]
    if(centrality=="mean") Centrality.feature <- apply(CURVES[Sample, ], 2, sum)/lg                
    if(centrality=="median") Centrality.feature <- median.npfda(SEMIMETRIC[Sample, Sample])
    MM0 <- rbind(Modal.curve, Centrality.feature, rep(0, p))
    SM.MM0 <- sm(MM0, MM0, ...)
    Hetero.index[e] <- SM.MM0[1,2]/(SM.MM0[1,3]+SM.MM0[2,3])
  }
  return(sum(Hetero.index)/nss)
}
#####################################################################
#####################################################################
#####################################################################
classif.split <- function(PROBCURVES, index, Group)
{
  #################################################################
  # Returns a partition of "Group"
  #   "PROBCURVES" contains the probability curves (row by row)
  #   "index" is the index correponding to the optimal bandwidth
  #           (see function "classif.bw")
  #   "Group" gives the unit numbers of the concerned group
  # Returns a list allowing of building the partition:
  #   if the density function of the probability points admits at 
  #   least one local minimum, then it is possible to split "Group" 
  #   into subgroups contained in the returned list.
  #################################################################
  Prob.points <- PROBCURVES[Group, index]
  pp.start <- min(Prob.points)
  pp.end <- max(Prob.points)
  subgroup <- list()
  if(pp.end > pp.start){
    est <- density(Prob.points, bw="bcv", from=pp.start, to=pp.end)
    Prob.lim <- est$x[rank.minima(est$y)]
    if(length(Prob.lim)==0){ 
      warning("Zero local minimum ==> this group is a terminal leaf of the classification tree")
    } else {
      kmax <- length(Prob.lim)
      subgroup[[1]] <- Group[Prob.points>Prob.lim[kmax]]
      subgroup[[kmax+1]] <- Group[Prob.points<= Prob.lim[1]]
      if(kmax>1){
        for(k in 1:(kmax-1))
          subgroup[[k+1]] <- Group[(Prob.points>Prob.lim[k]) & (Prob.points<=Prob.lim[k+1])]
      }
    }
  } else {
    warning("Zero local minimum ==> this group is a terminal leaf of the classification tree")
  }
  return(subgroup)
}
#####################################################################
#####################################################################
#####################################################################
entropy <- function(f, x)
{
  ##############################################################
  # Returns the entropy of the strictly positive function "f" 
  # known at points "x":
  # "f" is a numeric vector contained the values of the function
  # "x" is a numeric vector contained the values at which the 
  #     function is valued
  ##############################################################
  Diffx <- diff(x)
  f <- abs(f)
  flogf <- rep(0, length(f))
  logic <- (f!=0)
  flogf[logic] <- f[logic] * log(f[logic])
  n <- length(f)
  return(sum((flogf[-1] + flogf[ - n]) * Diffx)/2)
}
##############################################################################
##############################################################################
##############################################################################
##############################################################################
ffunopare.knn <- function(RESPONSES, CURVES, PRED, neighbour,..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional nonparametric prediction (regression) when both  
  # the response and the predictor are random curves via the functional kernel estimator. 
  # A global kNN bandwidth is given by the user.
  #    "RESPONSES" matrix containing the response curves (row by row) 
  #    "CURVES" matrix containing the explanatory curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "neighbour" number of nearest neighbours used by the estimator
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift" and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Bandwidths" vector containing the global data-driven bandwidths  
  #                 for each curve in the matrix "CURVES"
  #    "Cv" cross-validation criterion computed with "CURVES" and "RESPONSES"
  ################################################################
  if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
  SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  p1 <- ncol(SEMIMETRIC1)
  n1 <- ncol(SEMIMETRIC1)
  if(neighbour >= n1)
    stop(paste("try a smaller number of neighbour \n than ", neighbour))
  bandwidth.knn1 <- 0
  for(j in 1:p1) {
    Sem <- SEMIMETRIC1[, j]
    knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
    bandwidth.knn1[j] <- 0.5 * sum(knn.to.band)
  }
  KERNEL1 <- kernel(t(t(SEMIMETRIC1)/bandwidth.knn1))
  KERNEL1[KERNEL1 < 0] <- 0
  KERNEL1[KERNEL1 > 1] <- 0
  diag(KERNEL1) <- 0
  Denom1 <- apply(KERNEL1, 2, sum)
  NUM1 <- t(KERNEL1) %*% RESPONSES
  RESPONSES.estimated <- NUM1/Denom1
  Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p1)
  if(twodatasets) {
    SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    Bandwidth2 <- 0
    p2 <- ncol(SEMIMETRIC2)
    bandwidth.knn2 <- 0
    for(j in 1:p2) {
      Sem <- SEMIMETRIC2[, j]
      knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)
                                    ]]
      bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
    }
    KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
    KERNEL2[KERNEL2 < 0] <- 0
    KERNEL2[KERNEL2 > 1] <- 0
    Denom2 <- apply(KERNEL2, 2, sum)
    NUM2 <- t(KERNEL2) %*% RESPONSES
    RESPONSES.predicted <- NUM2/Denom2
    return(list(Estimated.values = RESPONSES.estimated, 
                Predicted.values = RESPONSES.predicted, Cv = 
                  Cv.estimated))
  }else {
    return(list(Estimated.values = RESPONSES.estimated, Cv = 
                  Cv.estimated))
  }
}
##############################################################################
##############################################################################
##############################################################################
##############################################################################
ffunopare.knn.gcv <- function(RESPONSES, CURVES, PRED, ..., Knearest=NULL, kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional nonparametric prediction (regression) when both  
  # the response and the predictor are random curves via the functional kernel estimator. 
  # A global bandwidth (i.e. number of neighbours) is selected by a 
  # cross-validation procedure.
  #    "RESPONSES" matrix containing the response curves (row by row) 
  #    "CURVES" matrix containing the explanatory curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "Knearest" vector giving the grid of nearest neighbours used in the 
  #               estimator ; if "Knearest"=NULL (default), a grid is 
  #               automatically built
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift" and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Bandwidths" vector containing the global data-driven bandwidths  
  #                 for each curve in the matrix "CURVES"
  #    "knearest.opt" optimal number of neighbours
  #    "Cv" cross-validation criterion computed with "CURVES" and "RESPONSES"
  ################################################################
  if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
  SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  if(is.null(Knearest)){
    step <- ceiling(n1/100)
    if(step == 0) step <- 1
    Knearest <- seq(from = 5, to = n1 %/% 2, by = step)	
    # the vector Knearest contains the sequence of the 
    # k-nearest neighbours used for computing the optimal bandwidth
  }
  kmax <- max(Knearest)	
  p <- ncol(CURVES)
  RESPONSES.estimated <- matrix(0, nrow = n1, ncol = p)
  Bandwidth.opt <- 0
  HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
  BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
  lKnearest <- length(Knearest)
  HAT.RESP <- array(0,c(n1,lKnearest,p))
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]	
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)	
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]	
    # BANDWIDTH[i, l-1] contains (dq(X_{j_l},X_i) + 
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/BANDWIDTH[i,  ]
    KNUM <- kernel(UMAT)
    KNUM[col(KNUM) > row(KNUM)] <- 0
    Kdenom <- apply(KNUM[Knearest,  ], 1, sum)
    WEIGHTS <- KNUM[Knearest,  ]/Kdenom
    Ind.curves <- Norm.order[2:(kmax + 1)]
    HAT.RESP[i,,] <- WEIGHTS %*% RESPONSES[Ind.curves,]
  }
  CRITARR <- array(0,c(n1,p,lKnearest))
  for(i in 1:n1){
    CRITARR[i,,] <- (t(HAT.RESP[i,,]) - RESPONSES[i,])^2
  }
  Criterium <- apply(CRITARR, 3, sum)
  index.opt <- order(Criterium)[1]
  RESPONSES.estimated <- HAT.RESP[, index.opt,]
  knearest.opt <- Knearest[index.opt]
  Bandwidth.opt <- BANDWIDTH[, knearest.opt]
  Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p)
  if(twodatasets) {
    SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(KERNEL, 2, sum)
    NUM <- t(KERNEL) %*% RESPONSES
    RESPONSES.predicted <- NUM/Denom
    return(list(Estimated.values = RESPONSES.estimated, 
                Predicted.values = RESPONSES.predicted, Bandwidths = 
                  Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                  Cv.estimated))
  }else {
    return(list(Estimated.values = RESPONSES.estimated, Bandwidths
                = Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                  Cv.estimated))
  }
}
##############################################################################
##############################################################################
##############################################################################
##############################################################################
ffunopare.knn.single <- function(RESPONSES, CURVES, PRED, neighbour,..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional nonparametric prediction (regression) when both  
  # the response and the predictor are random curves via the functional kernel estimator. 
  # A global kNN bandwidth is given by the user
  #    "RESPONSES" matrix containing the response curves (row by row) 
  #    "CURVES" matrix containing the explanatory curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "neighbour" number of nearest neighbours used by the estimator
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift" and "pca"
  # Returns a list containing:
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #    "knn" number of neighbours
  ################################################################
  if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
  SEMIMETRIC2 <- sm(CURVES, PRED, ...)
  kernel <- get(kind.of.kernel)
  p2 <- ncol(SEMIMETRIC2)
  bandwidth.knn2 <- 0
  for(j in 1:p2) {
    Sem <- SEMIMETRIC2[, j]
    knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
    bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
  }
  KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
  KERNEL2[KERNEL2 < 0] <- 0
  KERNEL2[KERNEL2 > 1] <- 0
  Denom2 <- apply(KERNEL2, 2, sum)
  NUM2 <- t(KERNEL2) %*% RESPONSES
  RESPONSES.predicted <- NUM2/Denom2
  return(list(Predicted.values = RESPONSES.predicted, knn=neighbour))
}
#####################################################################
#####################################################################
#####################################################################
fourier <- function(x, nbasis = n, period = span, nderiv = 0)
{
  #  Performed by J.O. Ramsay and available on its 
  #  website http://ego.psych.mcgill.ca/misc/fda which contains a lot 
  #  of functions for analyzing functional data in a different way:
  #  *Computes the NDERIV derivative of the Fourier series basis
  #   for NBASIS functions with period PERIOD, these being evaluated
  #   at values in vector X
  #  *Returns an N by NBASIS matrix of function values
  #  *Note:  The number of basis functions always odd.  If the argument
  #   NBASIS is even, it is increased by one.
  x <- as.vector(x)
  n <- length(x)
  onen <- rep(1, n)
  xrange <- range(x)
  span <- xrange[2] - xrange[1]
  if(nbasis <= 0)
    stop("NBASIS not positive")
  if(period <= 0)
    stop("PERIOD not positive")
  if(nderiv < 0)
    stop("NDERIV negative")
  if(2 * (nbasis %/% 2) == nbasis)
    nbasis <- nbasis + 1
  basis <- matrix(0, n, nbasis)
  omega <- (2 * pi)/period
  omegax <- omega * x
  if(nderiv == 0) {
    #  The fourier series itself is required.
    basis[, 1] <- 0.7071068
    j <- seq(2, nbasis - 1, 2)
    k <- j/2
    args <- outer(omegax, k)
    basis[, j] <- sin(args)
    basis[, j + 1] <- cos(args)
  }else {
    #  A derivative of the fourier series is required.
    basis[, 1] <- 0
    if(nderiv == floor(nderiv/2) * 2) {
      mval <- nderiv/2
      ncase <- 1
    }else {
      mval <- (nderiv - 1)/2
      ncase <- 2
    }
    j <- seq(2, nbasis - 1, 2)
    k <- j/2
    fac <- outer(onen, ((-1)^mval) * (k * omega)^nderiv)
    args <- outer(omegax, k)
    if(ncase == 1) {
      basis[, j] <- fac * sin(args)
      basis[, j + 1] <- fac * cos(args)
    }else {
      basis[, j] <- fac * cos(args)
      basis[, j + 1] <-  - fac * sin(args)
    }
  }
  basis <- basis/sqrt(period/2)
  return(basis)
}
#####################################################################
#####################################################################
#####################################################################
funopa.mode <- function(band, SEMIMETRIC, kind.of.kernel = "quadratic")
{
  #######################################################
  # Returns the rank of the modal curve in the curves 
  # sample used for computing the matrix "SEMIMETRIC":
  # "band" is the bandwidth used for estimating 
  #             the distribution of the curves sample
  # "SEMIMETRIC[i,j]" contains "d(Xi,Xj)" where "d(.,.)"
  #             is a fixed semimetric, Xi (resp. Xj) the 
  #             ith (resp. jth) curve of the curves sample
  ########################################################
  kernel <- get(kind.of.kernel)
  n <- nrow(SEMIMETRIC)
  diag(SEMIMETRIC) <- 0
  KERNEL <- kernel(SEMIMETRIC/band)
  KERNEL[KERNEL < 0] <- 0
  diag(KERNEL) <- 0
  density.estimate <- apply(KERNEL, 1, sum)/n
  rank.mode <- order(density.estimate)[n]
  return(rank.mode)
}
#####################################################################
#####################################################################
#####################################################################
funopadi.knn.lcv <- function(ind_train,data, kind.of.kernel = "quadratic", pairwise_dist,true_class)
{
  Classes<- true_class[ind_train]
  CURVES<- data[ind_train,]
  PRED<- data[-ind_train,]
  class_test<- true_class[-ind_train]
  ################################################################
  # Performs functional discrimination of a sample of curves when 
  # a categorical response is observed (supervised classification). 
  # A local bandwidth (i.e. local number of neighbours) is selected 
  # by a cross-validation procedure.
  #    "Classes" vector containing the categorical responses
  #              giving the group number for each curve in 
  #              the matrix CURVES (if nbclass is the number of 
  #              groups, "Classes" contains numbers 1,2,...,nbclass)
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.classnumber" vector containing estimated class membership 
  #                            for each curve of "CURVES"
  #    "Predicted.classnumber" if PRED different from CURVES, this vector 
  #                            contains predicted class membership for each 
  #                            curve of PRED
  #    "Bandwidths" vector containing the local data-driven bandwidths
  #                 for each curve in the matrix "CURVES"
  #    "Misclas" misclassification rate computed from estimated values and 
  #              observed values
  ################################################################
  Classes <- as.vector(Classes)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  SEMIMETRIC1 <- pairwise_dist[ind_train,ind_train]
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  step <- ceiling(n1/100)
  if(step == 0)
    step <- 1
  Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  kmax <- max(Knearest)
  # the vector Knearest contains the sequence of the 
  # k-nearest neighbours used for computing the optimal bandwidth
  Classes.estimated <- 0
  Bandwidth.opt <- 0
  nbclass <- max(Classes)
  BINARY <- matrix(0, n1, nbclass)
  for(g in 1:nbclass)
    BINARY[, g] <- as.numeric(Classes == g)
  HAT.PROB <- matrix(0, nrow = nbclass, ncol = length(Knearest))
  Knn1 <- 0
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]
    # Bandwidth[l-1] contains (dq(X_{j_l},X_i) + 
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    Bandwidth <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/Bandwidth
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0
    Ind.curves <- Norm.order[2:(kmax + 1)]
    for(g in 1:nbclass) {
      Ind.resp <- BINARY[Ind.curves, g]
      YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow
                     = T)
      HAT.PROB[g,  ] <- apply(YMAT[Knearest,  ] * KMAT[
        Knearest,  ], 1, sum)
    }
    Kmatsumbyrow <- apply(KMAT[Knearest,  ], 1, sum)
    HAT.PROB <- HAT.PROB/matrix(Kmatsumbyrow,nrow(HAT.PROB), ncol(HAT.PROB), byrow=T)
    Criterium <- t(rep(1, nbclass)) %*% (HAT.PROB - BINARY[i,  ])^2
    index <- order(as.vector(Criterium))[1]
    Knn1[i] <- Knearest[index]
    Classes.estimated[i] <- order(HAT.PROB[, index])[nbclass]
    Bandwidth.opt[i] <- Bandwidth[index]
  }
  Misclas.estimated <- sum(Classes.estimated != Classes)/n1
  if(twodatasets) {
    SEMIMETRIC2 <- pairwise_dist[ind_train,-ind_train]
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Sm2k.ord <- order(SEMIMETRIC2[, k])
      knn <- Knn1[Sm2k.ord[1]]
      Bandwidth2[k] <- sum(sort(Sm2k)[knn:(knn+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(as.matrix(KERNEL), 2, sum)
    PROB.PREDICTED <- matrix(0, nrow = n2, ncol = nbclass)
    for(g in 1:nbclass) {
      PROBKERNEL <- KERNEL * BINARY[, g]
      PROB.PREDICTED[, g] <- apply(as.matrix(PROBKERNEL), 2, 
                                   sum)/Denom
    }
    Classes.predicted <- as.vector((PROB.PREDICTED == apply(
      PROB.PREDICTED, 1, max)) %*% (1:nbclass))
    Misclas<- sum(Classes.predicted != class_test)/n2
    return(Misclas.estimated)
    #return(list(Estimated.classnumber = Classes.estimated, 
                #Predicted.classnumber = Classes.predicted, Bandwidths
                #= Bandwidth.opt, Misclas = Misclas.estimated))
    
  }else {
    return(list(Estimated.classnumber = Classes.estimated, 
                Bandwidths = Bandwidth.opt, Misclas = Misclas.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.kernel <- function(Response, CURVES, PRED, bandwidth, ..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A global bandwidth is considered  without automatic selection. 
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "bandwidth" the value of the bandwidth
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "band" value of the current bandwidth
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  KERNEL1 <- kernel(SEMIMETRIC1/bandwidth)
  KERNEL1[KERNEL1 < 0] <- 0
  KERNEL1[KERNEL1 > 1] <- 0
  diag(KERNEL1) <- 0
  RESPKERNEL1 <- KERNEL1 * Response
  Denom1 <- apply(KERNEL1, 2, sum)
  if(sum(Denom1 == 0) > 0)
    stop(paste("try a greater bandwidth than ", bandwidth))
  Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
  Mse.estimated <- sum((Response.estimated - Response)^2)/length(Response
  )
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    KERNEL2 <- kernel(SEMIMETRIC2/bandwidth)
    KERNEL2[KERNEL2 < 0] <- 0
    KERNEL2[KERNEL2 > 1] <- 0
    Denom2 <- apply(KERNEL2, 2, sum)
    if(sum(Denom2 == 0) > 0)
      stop(paste(
        "try a greater bandwidth \n                                                       than ",
        bandwidth))
    RESPKERNEL2 <- KERNEL2 * Response
    Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, band = bandwidth,
                Mse = Mse.estimated))
  }else {
    return(list(Estimated.values = Response.estimated, band = 
                  bandwidth, Mse = Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.kernel.cv <- function(Response, CURVES, PRED, ..., kind.of.kernel = "quadratic", semimetric = "deriv", h.range = NULL)
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A global bandwidth is automatically selected with a 
  # cross-validation procedure. 
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  #    "h.range" a vector of length 2 giving the range for the bandwidth. 
  #              By default, the procedure defines a sequence of candidates
  #              of bandwidths according to the values of the matrix CURVES 
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "hopt" value of the optimal bandwidth
  #    "hseq" the used sequence of possible bandwidths
  #    "Mse" mean squared error between estimated values 
  #          and observed values 
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  Semimetric1 <- SEMIMETRIC1[row(SEMIMETRIC1) > col(SEMIMETRIC1)]
  if(is.null(h.range)) {
    h.seq <- quantile(Semimetric1, seq(0.05, 0.5, length = 20))
  }else {
    h.seq <- seq(h.range[1], h.range[2], length = 20)
  }
  kernel <- get(kind.of.kernel)
  count1 <- 0
  count2 <- 0
  Mse <- 0
  h.seq.corrected <- 0
  h.seq.length <- length(h.seq)
  h <- h.seq[1]
  while(count1 < h.seq.length) {
    count1 <- count1 + 1
    h <- 1.1 * h
    KERNEL1 <- kernel(SEMIMETRIC1/h)
    diag(KERNEL1) <- 0
    KERNEL1[KERNEL1 < 0] <- 0
    KERNEL1[KERNEL1 > 1] <- 0
    Denom1 <- apply(KERNEL1, 2, sum)
    Logic <- (Denom1 == 0)
    while(sum(Logic) >= 1) {
      h <- 1.1 * h
      KERNEL1 <- kernel(SEMIMETRIC1/h)
      diag(KERNEL1) <- 0
      KERNEL1[KERNEL1 < 0] <- 0
      KERNEL1[KERNEL1 > 1] <- 0
      Denom1 <- apply(KERNEL1, 2, sum)
      Logic <- (Denom1 == 0)
    }
    RESPKERNEL1 <- KERNEL1 * Response
    Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
    count2 <- count2 + 1
    h.seq.corrected[count2] <- h
    Mse[count2] <- sum((Response.estimated - Response)^2)
  }
  index.opt <- order(Mse)[1]
  h.opt <- h.seq.corrected[index.opt]
  KERNEL1 <- kernel(SEMIMETRIC1/h.opt)
  KERNEL1[KERNEL1 < 0] <- 0
  KERNEL1[KERNEL1 > 1] <- 0
  diag(KERNEL1) <- 0
  Denom1 <- apply(KERNEL1, 2, sum)
  RESPKERNEL1 <- KERNEL1 * Response
  Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
  Mse.estimated <- sum((Response.estimated - Response)^2)/length(Response)
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    KERNEL2 <- kernel(SEMIMETRIC2/h.opt)
    KERNEL2[KERNEL2 < 0] <- 0
    KERNEL2[KERNEL2 > 1] <- 0
    Denom2 <- apply(KERNEL2, 2, sum)
    Logic <- (Denom2 == 0)
    while(sum(Logic) >= 1) {
      h.opt <- 1.1 * h.opt
      KERNEL2 <- kernel(SEMIMETRIC2/h.opt)
      diag(KERNEL2) <- 0
      KERNEL2[KERNEL2 < 0] <- 0
      KERNEL2[KERNEL2 > 1] <- 0
      Denom2 <- apply(KERNEL2, 2, sum)
      Logic <- (Denom2 == 0)
    }
    RESPKERNEL2 <- KERNEL2 * Response
    Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, hopt = h.opt, 
                hseq = h.seq.corrected, Mse = Mse.estimated))
  }else {
    return(list(Estimated.values = Response.estimated, hopt = h.opt,
                hseq = h.seq.corrected, Mse = Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.knn <- function(Response, CURVES, PRED, neighbour, ..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A bandwidth corresponding to number of neighbours has to be given.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "neighbour" number of neighbours fixed for computing the
  #                functional kernel estimator.
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "kNN" value of the current argument "neighbour".
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  p1 <- ncol(SEMIMETRIC1)
  n1 <- nrow(SEMIMETRIC1)
  if(neighbour >= n1)
    stop(paste("try a smaller number of neighbour \n                               than ", neighbour))
  bandwidth.knn1 <- 0
  for(j in 1:p1) {
    Sem <- SEMIMETRIC1[, j]
    knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
    bandwidth.knn1[j] <- 0.5 * sum(knn.to.band)
  }
  KERNEL1 <- kernel(t(t(SEMIMETRIC1)/bandwidth.knn1))
  KERNEL1[KERNEL1 < 0] <- 0
  KERNEL1[KERNEL1 > 1] <- 0
  diag(KERNEL1) <- 0
  RESPKERNEL1 <- KERNEL1 * Response
  Denom1 <- apply(KERNEL1, 2, sum)
  Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    p2 <- ncol(SEMIMETRIC2)
    bandwidth.knn2 <- 0
    for(j in 1:p2) {
      Sem <- SEMIMETRIC2[, j]
      knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)
                                    ]]
      bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
    }
    KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
    KERNEL2[KERNEL2 < 0] <- 0
    KERNEL2[KERNEL2 > 1] <- 0
    Denom2 <- apply(KERNEL2, 2, sum)
    RESPKERNEL2 <- KERNEL2 * Response
    Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, knn = neighbour, 
                Mse = Mse.estimated))
  }
  else {
    return(list(Estimated.values = Response.estimated, knn = 
                  neighbour, Mse = Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.knn.gcv <- function(Response, CURVES, PRED, ..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A global bandwidth (i.e. a number of neighbours) is selected by a 
  # cross-validation procedure.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Bandwidths" vector containing the global data-driven bandwidths  
  #                 for each curve in the matrix "CURVES"
  #    "knearest.opt" optimal number of neighbours
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  step <- ceiling(n1/100)
  if(step == 0)
    step <- 1
  Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  kmax <- max(Knearest)	
  # the vector Knearest contains the sequence of the 
  # k-nearest neighbours used for computing the optimal bandwidth
  Response.estimated <- 0
  Bandwidth.opt <- 0
  HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
  BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]	
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)	
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]	
    # BANDWIDTH[i, l-1] contains (dq(X_{j_l},X_i) + 
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/BANDWIDTH[i,  ]
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0
    Ind.curves <- Norm.order[2:(kmax + 1)]
    Ind.resp <- Response[Ind.curves]
    YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
    HAT.RESP[i,  ] <- apply(YMAT[Knearest,  ] * KMAT[Knearest,  ], 
                            1, sum)/apply(KMAT[Knearest,  ], 1, sum)
  }
  CRITERIUM <- (HAT.RESP - Response)^2
  Criterium <- apply(CRITERIUM, 2, sum)
  index.opt <- order(Criterium)[1]
  Response.estimated <- HAT.RESP[, index.opt]
  knearest.opt <- Knearest[index.opt]
  Bandwidth.opt <- BANDWIDTH[, knearest.opt]
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(KERNEL, 2, sum)
    RESPKERNEL <- KERNEL * Response
    Response.predicted <- apply(RESPKERNEL, 2, sum)/Denom
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, Bandwidths = 
                  Bandwidth.opt, knearest.opt = knearest.opt, Mse = 
                  Mse.estimated))
  }else {
    return(list(Estimated.values = Response.estimated, Bandwidths
                = Bandwidth.opt, knearest.opt = knearest.opt, Mse = 
                  Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.knn.lcv <- function(Response, CURVES, PRED, ..., kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional prediction (regression) of a scalar response 
  # from a sample of curves via the functional kernel estimator. 
  # A local bandwidth (i.e. local number of neighbours) is selected 
  # by a cross-validation procedure.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Bandwidths" vector containing the local data-driven bandwidths
  #                 for each curve in the matrix "CURVES"
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  if(semimetric == "mplsr")
    SEMIMETRIC1 <- sm(Response, CURVES, CURVES, ...)
  else SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
  kernel <- get(kind.of.kernel)
  n1 <- ncol(SEMIMETRIC1)
  step <- ceiling(n1/100)
  if(step == 0)
    step <- 1
  Knearest <- seq(from = 10, to = n1 %/% 2, by = step)
  kmax <- max(Knearest)
  # the vector Knearest contains the sequence of the 
  # k-nearest neighbours used for computing the optimal bandwidth
  Response.estimated <- 0
  Bandwidth.opt <- 0
  Knn1 <- 0
  for(i in 1:n1) {
    Norm.diff <- SEMIMETRIC1[, i]
    # "norm.order" gives the sequence k_1, k_2,... such that
    # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
    Norm.order <- order(Norm.diff)
    # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
    # dq(X_{j_{kamx+2}},X_i)
    zz <- sort(Norm.diff)[2:(kmax + 2)]
    # Bandwidth[l-1] contains (dq(X_{j_l},X_i) + 
    # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
    Bandwidth <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
    z <- zz[ - (kmax + 1)]
    ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
    UMAT <- ZMAT/Bandwidth
    KMAT <- kernel(UMAT)
    KMAT[col(KMAT) > row(KMAT)] <- 0
    Ind.curves <- Norm.order[2:(kmax + 1)]
    Ind.resp <- Response[Ind.curves]
    YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
    Hat.resp <- apply(YMAT[Knearest,  ] * KMAT[Knearest,  ], 1, sum
    )/apply(KMAT[Knearest,  ], 1, sum)
    Criterium <- abs(Hat.resp - Response[i])
    index <- order(Criterium)[1]
    Knn1[i] <- Knearest[index]
    Response.estimated[i] <- Hat.resp[index]
    Bandwidth.opt[i] <- Bandwidth[index]
  }
  Mse.estimated <- sum((Response.estimated - Response)^2)/n1
  if(twodatasets) {
    if(semimetric == "mplsr")
      SEMIMETRIC2 <- sm(Response, CURVES, PRED, ...)
    else SEMIMETRIC2 <- sm(CURVES, PRED, ...)
    Bandwidth2 <- 0
    n2 <- ncol(SEMIMETRIC2)
    for(k in 1:n2) {
      Sm2k <- SEMIMETRIC2[, k]
      Sm2k.ord <- order(SEMIMETRIC2[, k])
      knn <- Knn1[Sm2k.ord[1]]
      Bandwidth2[k] <- sum(sort(Sm2k)[knn:(knn+1)])*0.5
    }
    KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
    KERNEL[KERNEL < 0] <- 0
    KERNEL[KERNEL > 1] <- 0
    Denom <- apply(as.matrix(KERNEL), 2, sum)
    RESPKERNEL <- KERNEL * Response
    Response.predicted <- apply(as.matrix(RESPKERNEL), 2, sum)/
      Denom
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, Bandwidths = 
                  Bandwidth.opt, Mse = Mse.estimated))
  }
  else {
    return(list(Estimated.values = Response.estimated, Bandwidths
                = Bandwidth.opt, Mse = Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.mode.lcv <- function(Response, CURVES, PRED, ..., Knearest = NULL, kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional prediction of a scalar response from a 
  # sample of curves by computing the functional conditional mode. 
  # A local bandwidth (i.e. local number of neighbours) is selected 
  # by a ``trivial'' cross-validation procedure.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "Knearest"  vector giving the the sequence of successive authorized 
  #                integers for the smoothing parameters. By default 
  #                (i.e. Knearest=NULL), the vector Knearest contains a 
  #                sequence of 10 integers taking into account card(I1).
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" vector containing estimated reponses for 
  #                        each curve of "CURVES"
  #    "Predicted.values" if PRED different from CURVES, this vector 
  #                       contains predicted responses for each 
  #                       curve of PRED
  #    "Response.values"  vector of length card(I2) such that, for all 
  #                       i in the set I2, Response.values[i]=yi 
  #                       (i.e. observed responses corresponding to the 
  #                       second learning subsample).
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  kernel <- get(kind.of.kernel)
  llearn <- nrow(CURVES)
  Learn1 <- seq(2, llearn, by=2)
  llearn1 <- length(Learn1)
  LEARN1 <- CURVES[Learn1,  ]
  LEARN2 <- CURVES[ - Learn1,  ]
  if(semimetric == "mplsr"){
    SMLEARN1 <- sm(Response, LEARN1, LEARN1, ...)
    SMLEARN12 <- sm(Response, LEARN1, LEARN2, ...)
  } else {
    SMLEARN1 <- sm(LEARN1, LEARN1, ...)
    SMLEARN12 <- sm(LEARN1, LEARN2, ...)
  }
  SML12.SOR <- apply(SMLEARN12, 2, sort)
  Resp1 <- Response[Learn1]
  Resp2 <- Response[ - Learn1]
  Resp.range <- range(Response)
  Response.grid <- seq(from = Resp.range[1] * 0.9, to = Resp.range[2] * 
                         1.1, length = 100)	
  # RESPMETRIC[i,j]=yi-yj with i in Response.grid and j in LEARN1 
  RESPMETRIC <- outer(Response.grid, Resp1, "-")
  RESPMET.SOR <- t(apply(abs(RESPMETRIC), 1, sort))
  llearn2 <- nrow(LEARN2)
  lgrid <- length(Response.grid)
  if(is.null(Knearest)) {
    Knearest.min <- max(ceiling(llearn1 * 0.05), 10)
    Knearest.max <- ceiling(llearn1 * 0.25)
    if(Knearest.max <= Knearest.min){
      Knearest.min <- ceiling(llearn1 * 0.05)
    }
    step <- ceiling((Knearest.max - Knearest.min)/10)
    Knearest <- seq(Knearest.min, Knearest.max, by = step)
  }
  lknearest <- length(Knearest)
  BANDL12.CUR <- 0.5 * (SML12.SOR[Knearest,  ] + SML12.SOR[Knearest + 1,  
                                                           ])
  BAND.RESP <- 0.5 * (RESPMET.SOR[, Knearest] + RESPMET.SOR[, Knearest + 
                                                              1])
  CV <- matrix(0, nrow = lknearest^2, ncol = llearn2)
  MODE <- matrix(0, nrow = lknearest^2, ncol = llearn2)
  count1 <- 0
  count2 <- 0
  for(kk in Knearest) {
    count2 <- count2 + 1
    ARG <- t(t(SMLEARN12)/BANDL12.CUR[count2,  ])
    KERNEL.CURVES <- kernel(ARG)
    KERNEL.CURVES[KERNEL.CURVES < 0] <- 0
    KERNEL.CURVES[KERNEL.CURVES > 1] <- 0
    Denom <- apply(KERNEL.CURVES, 2, sum)
    count3 <- 0
    for(hh in Knearest) {
      count1 <- count1 + 1
      count3 <- count3 + 1
      KERNEL.RESP <- apply(abs(RESPMETRIC)/BAND.RESP[, count3
                                                     ], 1, kernel)
      KERNEL.RESP[KERNEL.RESP < 0] <- 0
      KERNEL.RESP[KERNEL.RESP > 1] <- 0
      DENSITY.ESTIMATE <- (t(KERNEL.CURVES)/Denom) %*% (
        KERNEL.RESP/BAND.RESP[, count3])
      Ind.mode <- apply(DENSITY.ESTIMATE, 1, order)[lgrid,  ]
      MODE[count1,  ] <- Response.grid[Ind.mode]
      CV[count1,  ] <- (Resp2 - MODE[count1,  ])^2
    }
  }
  Ind.knearest.opt <- apply(CV, 2, order)[1,  ]
  IND.OPT <- cbind(Ind.knearest.opt, 1:llearn2)
  Response.estimated <- MODE[IND.OPT]
  Mse.estimated <- sum(CV[IND.OPT])/llearn2
  if(twodatasets) {
    if(semimetric == "mplsr")
      SMLEARN2NEW <- sm(Response, LEARN2, PRED, ...)
    else SMLEARN2NEW <- sm(LEARN2, PRED, ...)
    Order.new <- apply(SMLEARN2NEW, 2, order)[1,  ]
    Response.predicted <- Response.estimated[Order.new]
    return(list(Estimated.values = Response.estimated, 
                Predicted.values = Response.predicted, Response.values
                = Resp2, Mse = Mse.estimated))
  }
  else {
    return(list(Estimated.values = Response.estimated, 
                Response.values = Resp2, Mse = Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
funopare.quantile.lcv <- function(Response, CURVES, PRED, ..., alpha = c(0.05, 0.5, 0.95), Knearest = NULL, kind.of.kernel = "quadratic", semimetric = "deriv")
{
  ################################################################
  # Performs functional prediction of a scalar response from a 
  # sample of curves by computing the functional conditional mode. 
  # A local bandwidth (i.e. local number of neighbours) is selected 
  # by a ``trivial'' cross-validation procedure.
  #    "Response" vector containing the observations of the scalar 
  #               response
  #    "CURVES" matrix containing the curves dataset (row by row) 
  #             used for the estimating stage
  #    "PRED" matrix containing new curves stored row by row
  #           used for computing predictions
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "alpha"  vector giving the quantiles to be computed. By default, 
  #             the 5-percentile, median and 95-percentile are computed.
  #    "Knearest"  vector giving the the sequence of successive authorized 
  #                integers for the smoothing parameters. By default 
  #                (i.e. Knearest=NULL), the vector Knearest contains a 
  #                sequence of 10 integers taking into account card(I1).
  #    "kind.of.kernel" the kernel function used for computing of 
  #                     the kernel estimator; you can choose 
  #                     "indicator", "triangle" or "quadratic (default)
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  # Returns a list containing:
  #    "Estimated.values" a  card(I2)-by-length(alpha) matrix whose 
  #                       columns gives the corresponding estimated 
  #                       conditional quantiles, for all curves in the 
  #                       second learning subsample I2, 
  #    "Predicted.values" if PRED different from CURVES, this matrix 
  #                       contains predicted conditional quantiles 
  #                       for each curve of PRED
  #    "Response.values"  vector of length the size of the second 
  #                       learning subsample giving the corresponding  
  #                       observed responses.
  #    "Mse" mean squared error between estimated values and 
  #          observed values
  ################################################################
  Response <- as.vector(Response)
  Logicalpha <- alpha==0.5
  nomedian <- !as.logical(sum(Logicalpha))
  if(nomedian) stop("Please, add in the argument alpha the median (i.e 0.5)")
  lalpha <- length(alpha)
  testalpha <- lalpha > 1
  if(testalpha) posmedian <- (1:lalpha)[Logicalpha] 
  if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
  onerow <- nrow(PRED) == 1
  testfordim <- sum(dim(CURVES)==dim(PRED))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  kernel <- get(kind.of.kernel)
  int.kernel <- get(paste("integrated.", kind.of.kernel, sep = ""))
  llearn <- nrow(CURVES)
  Learn1 <- seq(2, llearn, by=2)
  llearn1 <- length(Learn1)
  LEARN1 <- CURVES[Learn1,  ]
  LEARN2 <- CURVES[ - Learn1,  ]
  if(semimetric == "mplsr"){
    SMLEARN1 <- sm(Response, LEARN1, LEARN1, ...)
    SMLEARN12 <- sm(Response, LEARN1, LEARN2, ...)
  } else {
    SMLEARN1 <- sm(LEARN1, LEARN1, ...)
    SMLEARN12 <- sm(LEARN1, LEARN2, ...)
  }
  SML12.SOR <- apply(SMLEARN12, 2, sort)
  Resp1 <- Response[Learn1]
  Resp2 <- Response[ - Learn1]
  Resp.range <- range(Response)
  Response.grid <- seq(from = Resp.range[1] * 0.9, to = Resp.range[2] * 
                         1.1, length = 100)	
  # RESPMETRIC[i,j]=yi-yj with i in Response.grid and j in LEARN1 
  RESPMETRIC <- outer(Response.grid, Resp1, "-")
  RESPMET.SOR <- t(apply(abs(RESPMETRIC), 1, sort))
  llearn2 <- nrow(LEARN2)
  lgrid <- length(Response.grid)
  if(is.null(Knearest)) {
    Knearest.min <- max(ceiling(llearn1 * 0.05), 10)
    Knearest.max <- ceiling(llearn1 * 0.25)
    if(Knearest.max <= Knearest.min){
      Knearest.min <- ceiling(llearn1 * 0.05)
    }
    step <- ceiling((Knearest.max - Knearest.min)/10)
    Knearest <- seq(Knearest.min, Knearest.max, by = step)
  }
  lknearest <- length(Knearest)
  BANDL12.CUR <- 0.5 * (SML12.SOR[Knearest,  ] + SML12.SOR[Knearest + 1, ])
  BAND.RESP <- 0.5 * (RESPMET.SOR[, Knearest] + RESPMET.SOR[, Knearest + 
                                                              1])
  CV <- matrix(0, nrow = lknearest^2, ncol = llearn2)
  if(testalpha){
    QUANT <- array(0, dim = c(lknearest^2, llearn2, lalpha))
    
  } else {
    QUANT <- matrix(0, nrow = lknearest^2, ncol = llearn2)
  }
  count <- 0
  for(kk in 1:lknearest) {
    ARG <- t(t(SMLEARN12)/BANDL12.CUR[kk,])
    KERNEL.CURVES <- kernel(ARG)
    KERNEL.CURVES[KERNEL.CURVES < 0] <- 0
    KERNEL.CURVES[KERNEL.CURVES > 1] <- 0
    Denom <- apply(KERNEL.CURVES, 2, sum)
    for(hh in 1:lknearest) {
      count <- count + 1
      IKERNEL.RESP <- apply(RESPMETRIC/BAND.RESP[, hh], 1,
                            int.kernel)
      CDF.EST <- (t(KERNEL.CURVES)/Denom) %*% IKERNEL.RESP
      if(testalpha){
        for(ii in 1:lalpha){
          Ind.quant <- apply(CDF.EST < alpha[ii], 1, sum)
          Ind.quant[Ind.quant==0] <- 1
          QUANT[count,,ii] <- Response.grid[Ind.quant]
        }
        CV[count,] <- (Resp2 - QUANT[count,,posmedian])^2
      } else {
        Ind.quant <- apply(CDF.EST < alpha, 1, sum)
        Ind.quant[Ind.quant==0] <- 1
        QUANT[count, ] <- Response.grid[Ind.quant]
        CV[count, ] <- (Resp2 - QUANT[count,  ])^2
      }
    }
  }
  Ind.knearest.opt <- apply(CV, 2, order)[1,  ]
  IND.OPT <- cbind(Ind.knearest.opt, 1:llearn2)
  if(testalpha){
    Indkopt <- rep(Ind.knearest.opt, lalpha)
    Units <- rep(1:llearn2, lalpha) 
    Typeofest <- sort(rep(1:lalpha,llearn2))
    RESPONSE.ESTIMATED <- matrix(QUANT[cbind(Indkopt,Units,Typeofest)], nrow=llearn2, byrow=F)
    dimnames(RESPONSE.ESTIMATED) <- list(NULL, as.character(alpha))
  } else {
    RESPONSE.ESTIMATED <- QUANT[IND.OPT]
  }
  Mse.estimated <- sum(CV[IND.OPT])/llearn2
  if(twodatasets) {
    if(semimetric == "mplsr")
      SMLEARN2NEW <- sm(Response, LEARN2, PRED, ...)
    else SMLEARN2NEW <- sm(LEARN2, PRED, ...)
    Order.new <- apply(SMLEARN2NEW, 2, order)[1,  ]
    if(testalpha){
      if(onerow){
        RESPONSE.PREDICTED <- as.matrix(t(RESPONSE.ESTIMATED[Order.new,]))
      } else {
        RESPONSE.PREDICTED <- as.matrix(RESPONSE.ESTIMATED[Order.new,])
        dimnames(RESPONSE.PREDICTED) <- list(NULL, as.character(alpha))}
    } else {
      RESPONSE.PREDICTED <- RESPONSE.ESTIMATED[Order.new]
    }
    return(list(Estimated.values = RESPONSE.ESTIMATED, 
                Predicted.values = RESPONSE.PREDICTED, Response.values
                = Resp2, Mse = Mse.estimated))
  }
  else {
    return(list(Estimated.values = RESPONSE.ESTIMATED, 
                Response.values = Resp2, Mse = Mse.estimated))
  }
}
#####################################################################
#####################################################################
#####################################################################
kernel_smoother = function(bwsmooth, Std, CONTAMINATED_CURVES, Grid)
{
  ################################################################
  # Performs the denoising step on the contaminated profiles. A local 
  # bandwidth (proportional to the standard deviation derived from 
  # the measurement point where the smoother is valued) is used.   
  ################################################################
  #    "bwsmooth" given bandwidth
  #    "Std" vector containing standard deviations computed at the sampled points
  #    "CONTAMINATED_CURVES" matrix containing the contamined functional predictors (row by row)
  #    "Grid" vector of measurement points (the same for all curves)
  # 
  # Returns the relative k-fold cross-validation
  ################################################################
  GRID = outer(Grid, Grid, "-")
  ARGUMENT = GRID / (bwsmooth * Std)
  KERNEL = 1 - (ARGUMENT)^2
  KERNEL[ARGUMENT < -1] = 0
  KERNEL[ARGUMENT > 1] = 0
  Denominator = apply(KERNEL, 1, sum)
  WEIGHTS = KERNEL / Denominator
  OUTPUTS = CONTAMINATED_CURVES %*% WEIGHTS
  return(OUTPUTS)
}
#####################################################################
#####################################################################
#####################################################################
hshift <- function(x,y, grid)
{
  ####################################################################
  # Returns the "horizontal shifted proximity" between two discretized 
  # curves "x" and "y" (vectors of same length). 
  # The user has to choose a "grid".
  #####################################################################
  lgrid <- length(grid)
  a <- grid[1]
  b <- grid[lgrid]
  rang <- b - a
  lagmax <- floor(0.2 * rang)
  integrand <- (x-y)^2
  Dist1 <- sum(integrand[-1] + integrand[-lgrid])/(2 * rang)
  Dist2 <- Dist1
  for(i in 1:lagmax){
    xlag <- x[-(1:i)]
    xforward <- x[-((lgrid-i+1):lgrid)]
    ylag <- y[-(1:i)]
    yforward <- y[-((lgrid-i+1):lgrid)]
    integrand1 <- (xlag-yforward)^2
    integrand2 <- (xforward-ylag)^2
    lintegrand <- length(integrand1)
    rescaled.range <- 2 * (rang - 2 * i)
    Dist1[i+1] <- sum(integrand1[-1] + integrand1[-lintegrand])/rescaled.range
    Dist2[i+1] <- sum(integrand2[-1] + integrand2[-lintegrand])/rescaled.range
  }
  lag1 <- (0:lagmax)[order(Dist1)[1]]
  lag2 <- (0:lagmax)[order(Dist2)[1]]
  distmin1 <- min(Dist1)
  distmin2 <- min(Dist2)
  if(distmin1 < distmin2){
    distmin <- distmin1
    lagopt <- lag1
  }else{
    distmin <- distmin2
    lagopt <- lag2		
  }
  return(list(dist=sqrt(distmin),lag=lagopt))
}
#####################################################################
#####################################################################
#####################################################################
indicator <- function(u)
{
  Logic0 <- u<0
  Logic1 <- u>1
  Logic01 <- as.logical((1-Logic0) * (1-Logic1))
  u[Logic0] <- 0
  u[Logic1] <- 0
  u[Logic01] <- 1
  return(u)
}
#####################################################################
#####################################################################
#####################################################################
integrated.quadratic <- function(u)
{
  #  integrated quadratic kernel
  result <- u
  Logic0 <- (u <= -1)
  Logic1 <- (u >= 1)
  Logic2 <- (u > -1 & u < 1)
  result[Logic0] <- 0
  result[Logic1] <- 1
  Uval <- result[Logic2]
  result[Logic2] <- 0.75 * Uval * (1 - (Uval^2)/3) + 0.5
  return(result)
}
#####################################################################
#####################################################################
#####################################################################
integrated.triangle <- function(u)
{
  #  integrated triangle kernel
  result <- u
  Logic0 <- (u <= -1 | u >= 1)
  Logic1 <- (u > -1 & u < 0)
  Logic2 <- (u >=0 & u < 1)
  Uneg <- u[Logic1]
  Upos <- u[Logic2]
  result[Logic1] <-  Uneg + (Uneg^2)/2 + 0.5
  result[Logic2] <-  Upos - (Upos^2)/2 + 0.5
  result[Logic0] <- 0
  return(result)
}
#####################################################################
#####################################################################
############### IDENTIFYING INFLUENCE POINTS THROUGH ################
############### LOCAL LINEAR MULTIVARIATE REGRESSION ################
############### WITH SIMPLIFIED BANDWIDTHS:          ################
############### STEPWISE ALGORITHM                   ################
#####################################################################
#####################################################################
#####################################################################
mpdp <- function(CURVES, Response, nbmax=5, nbbw=5, Grid=(1:ncol(CURVES)), pcvpar=1, kind.of.kernel="quadratic2")
  ##############################################################################
# Select most-predictive design points ("mpdp") when regressing a saclar response 
# "Response" on functional data predictors "CURVES"
# "CURVES" matrix (n rows and p columns) containing the functional data predictors 
#          X1, X2,...,Xn (n curves) stored row by row: 
#          row 1: X1(t1),...,X1(tp) 
#          row 2: X2(t1),...,X2(tp) 
#          .......................
#          row p: Xn(t1),...,Xn(tp) 
#          it is worth noting that this sample of curves must be observed at the 
#          same p design points t1,...,tp
# "Response" vector of length n containing the corresponding responses Y1,...,Yn 
# "nbmax" maximum number of selected design points (default value: 5)
# "nbbw" size of the grid of bandwidths (h1,h2,..,hnbbw) used to choose the optimal 
#          one. From "nbbw", a grid of bandwidths is built in terms of k-nearest
#          neighbours: 
#          h1 takes into account 2% of the sample (for building the kernel estimator)
#          ......................................  
#          hk takes into account [2+{(50-2)/(nbbw-1)}*(k-1)]% of the sample
#          ......................................  
#          hnbbw takes into account 50% of the sample
#          (default value 5 ==> 2%, 14%, 26%, 38%, 50%)
# "Grid" vector containg a subset of design points in which one looks for the 
#          "nbmax" most-predicitve ones (default value: the whole set of 
#          design points is considered). 
# "pcvpar" parameter used in the penalized corssvalidation criterion (delta0 in 
#          the corresponding paper). The default value is 1. 
# "kind.of.kernel" gives the kernel used for computing the modes (default kernel 
#          is the quadratic one)          
##############################################################################
# Returns an object of class "mpdp": 
# "$Bwdselection"  vector identifying the most-predicitve design points retained 
#                    definitely after the backward deletion stage
# "$Fwdselection"  vector identifying the most-predicitve design points retained 
#                    just after the forward addition stage
# "$bandwidth"     optimal bandwidth (see paper for more details)
# "$Bwd.estimated.values" vector containing the leave-one-out local linear 
#                    estimations of the responses obtained by using the  
#                    selected most-predictive design points  
# "$Bwdpcv"        vector containing the values of the penalized crossvalidation 
#                    criterion at each step of the backward deletion stage 
# "$Fwdpcv"        vector containing the values of the penalized crossvalidation 
#                    criterion at each step of the forward addition stage 
#                          final penalized crossvalidation criterion
# "$fpcv.status"   "yes" means pcv reached its minimum at the forward addition stage;
#                  "no" means that the number of selected design points to reach the 
#                    minimum must be greater than "nbmax" (i.e. start again "mpdp" 
#                    with a larger "nbmax")
##############################################################################
{
  Selected.points <- NULL
  XX <- NULL
  kernel <- get(kind.of.kernel)
  nind <- nrow(CURVES)
  nind2 <- nind^2
  Ind11 <- rep(1:nind,nind)
  Ind12 <- rep(1:nind,rep(nind,nind))
  Gridnew <- Grid
  Gridpoints <- NULL
  poszero <- ((0:(nind-1))*nind)+(1:nind)
  Nbknn <- ceiling(seq(0.02,0.5,length=nbbw)*nind)
  varesp <- var(Response)
  Leaveoneout.estimated.values <- list()
  Meancurve <- apply(CURVES,2,mean)
  CCURVES <- t(t(CURVES)-Meancurve)
  Stdev <- sqrt(apply(CCURVES^2,2,mean))
  STDCURVES <- t(t(CCURVES)/Stdev)
  Bandwidths <- 0
  Fcv <- 0
  Fpcv <- 0
  minimum.reached.for.pcv <- "no"
  fpcv.ini <- sum((Response-mean(Response))^2)/nind
  ##########################################################
  # FORWARD ADDITION
  ##########################################################
  for(kk in 1:nbmax){
    cat(paste("Forward addition in progress"),sep="\n")
    if(length(Selected.points)!=0){Gridnew <- Gridnew[-selectedpoint]}
    lgridnew <- length(Gridnew)
    INVSUMOFSQUARES <- matrix(0,nrow=nbbw,ncol=lgridnew)
    ESTIMATIONS <- array(0,c(nbbw,lgridnew,nind))
    jj=0
    for(tt in Gridnew){
      jj=jj+1
      Xxtt <- STDCURVES[,tt]
      XXNEW <- cbind(XX, Xxtt)
      XXNEW.ALT.ROW <- XXNEW[Ind11,]
      XXNEW.REP.ROW <- XXNEW[Ind12,]
      D2 <- (XXNEW.ALT.ROW-XXNEW.REP.ROW)^2
      Prox <- sqrt(apply(as.matrix(D2),1,sum))
      PROX <- matrix(Prox, nrow=nind)
      Max <- apply(PROX,1,max) 
      diag(PROX)=Max 
      PROXSORTED <- apply(PROX,1,sort)
      Bw <- apply(PROXSORTED[Nbknn,],1,max)
      ll=0
      for(ll in 1:nbbw){
        Dist <- Prox/Bw[ll]
        Ker <- rep(0,nind2)
        Ker[-poszero] <- kernel(Dist[-poszero])
        KERNEL <- matrix(Ker, nrow=nind)
        Denom <- apply(KERNEL,1,sum)
        if(sum(Denom==0)>0){
          INVSUMOFSQUARES[ll,jj] <- 0
          next
        }else{
          NUM <- KERNEL %*% XXNEW
          XBAR <- NUM/Denom
          Ynum <- as.vector(KERNEL %*% Response) 
          Ybar <- Ynum/Denom
          T11 <- crossprod(t(KERNEL)*Response,XXNEW)
          T12 <- XBAR*Ynum
          T1 <- T11-T12
          for(ii in 1:nind){
            CXXNEW <- t(t(XXNEW)-XBAR[ii,])
            T2 <- crossprod(CXXNEW,(CXXNEW*KERNEL[ii,]))
            if(abs(det(T2))<1e-8){
              ESTIMATIONS[ll,jj,ii] <- Ybar[ii]
            }else{
              Bb <- solve(T2,T1[ii,])
              ESTIMATIONS[ll,jj,ii] <- Ybar[ii]+crossprod(Bb,CXXNEW[ii,])
            } 
          }
          INVSUMOFSQUARES[ll,jj] <- 1/sum((Response-ESTIMATIONS[ll,jj,])^2)
        }
      }
    }
    nr <- nrow(INVSUMOFSQUARES)
    Numcolofmax <- max.col(INVSUMOFSQUARES)
    numrowofmax <- which.max(INVSUMOFSQUARES[cbind(1:nr,Numcolofmax)])
    selectedpoint <- Numcolofmax[numrowofmax]
    Selected.points <- c(Selected.points,selectedpoint)
    Gridpoints <- c(Gridpoints,Gridnew[selectedpoint])
    Xxtt <- STDCURVES[,Gridnew[selectedpoint]]
    XX <- cbind(XX, Xxtt)
    XX.ALT.ROW <- XX[Ind11,]
    XX.REP.ROW <- XX[Ind12,]
    D2 <- (XX.ALT.ROW-XX.REP.ROW)^2
    Prox <- sqrt(apply(as.matrix(D2),1,sum))
    PROX <- matrix(Prox, nrow=nind)
    Max <- apply(PROX,1,max) 
    diag(PROX) <- Max 
    PROXSORTED <- apply(PROX,1,sort)
    Bw <- apply(PROXSORTED[Nbknn,],1,max)
    Bandwidths[kk] <- Bw[numrowofmax]
    Leaveoneout.estimated.values[[kk]] <- ESTIMATIONS[numrowofmax,Numcolofmax[numrowofmax],]
    Fcv[kk] <- sum((Response-Leaveoneout.estimated.values[[kk]])^2)/nind
    Fpcv[kk] <- Fcv[kk]*(1+kk*pcvpar/log(nind))
    if(kk==1) {
      if(Fpcv[kk]>fpcv.ini) {
        stop("No design point selected: stepwise algorithm over")
      } else {
        cat(paste("Design point",Gridnew[selectedpoint],"selected"),sep="\n")
      }
    } else {
      if(Fpcv[kk]>Fpcv[kk-1]) {
        minimum.reached.for.pcv <- "yes"
        break      
      } else {
        cat(paste("Design point",Gridnew[selectedpoint],"selected"),sep="\n")
      }
    }
  }
  if(minimum.reached.for.pcv=="yes") {
    Forwardselection <- Gridpoints[-kk]
    Fbwstd <- Bandwidths[kk-1]
    Fleaveoneout.estimated.values <- Leaveoneout.estimated.values[[kk-1]]
    Fpcv <- Fpcv[-kk]
  } else {
    cat("Minimum for PCV not reached; increase nbmax", sep="\n")
    Forwardselection <- Gridpoints
    Fbwstd <- Bandwidths[kk]
    Fleaveoneout.estimated.values <- Leaveoneout.estimated.values[[nbmax]]
  }
  ##########################################################
  # BACKWARD DELETION
  ##########################################################
  Selected.points <- NULL
  Bandwidths <- 0
  XX <- STDCURVES[,Forwardselection]
  Gridpoints <- NULL
  Gridnew <- Forwardselection
  lforwardselection <- length(Forwardselection)
  Bcv <- 0
  Bpcv <- Fpcv[lforwardselection]
  exitloop <- F
  for(kk in 1:(lforwardselection-1)){
    if(kk==1) {abbrev <- "st"}
    if(kk==2) {abbrev <- "nd"}
    if(kk==3) {abbrev <- "rd"}
    if(kk>3) {abbrev <- "th"}
    cat(paste("Backward deletion: is the deletion of a ",kk,abbrev," design point necessary?", sep=""),sep="\n")
    if(length(Selected.points)!=0){Gridnew <- Gridnew[-selectedpoint]}
    lgridnew <- length(Gridnew)
    SUMOFSQUARES <- matrix(0,nrow=lgridnew,ncol=nbbw)
    ESTIMATIONS <- array(0,c(nbbw,lgridnew,nind))
    jj=0
    for(tt in 1:lgridnew){
      jj=jj+1
      XXNEW <- XX[,-tt]
      XXNEW.ALT.ROW <- XXNEW[Ind11,]
      XXNEW.REP.ROW <- XXNEW[Ind12,]
      D2 <- (XXNEW.ALT.ROW-XXNEW.REP.ROW)^2
      Prox <- sqrt(apply(as.matrix(D2),1,sum))
      PROX <- matrix(Prox, nrow=nind)
      Max <- apply(PROX,1,max) 
      diag(PROX)=Max 
      PROXSORTED <- apply(PROX,1,sort)
      Bw <- apply(PROXSORTED[Nbknn,],1,max)
      ll=0
      for(ll in 1:nbbw){
        Dist <- Prox/Bw[ll]
        Ker <- rep(0,nind2)
        Ker[-poszero] <- kernel(Dist[-poszero])
        KERNEL <- matrix(Ker, nrow=nind)
        Denom <- apply(KERNEL,1,sum)
        if(sum(Denom==0)>0){
          SUMOFSQUARES[jj,ll] <- sum((Response-mean(Response))^2)
          next
        }else{
          NUM <- KERNEL %*% XXNEW
          XBAR <- NUM/Denom
          Ynum <- as.vector(KERNEL %*% Response) 
          Ybar <- Ynum/Denom
          T11 <- crossprod(t(KERNEL)*Response,XXNEW)
          T12 <- XBAR*Ynum
          T1 <- T11-T12
          for(ii in 1:nind){
            #cat(paste("ii ",ii),sep="\n")
            CXXNEW <- t(t(XXNEW)-XBAR[ii,])
            T2 <- crossprod(CXXNEW,(CXXNEW*KERNEL[ii,]))
            if(abs(det(T2))<1e-8){
              ESTIMATIONS[ll,jj,ii] <- Ybar[ii]
            }else{
              Bb <- solve(T2,T1[ii,])
              ESTIMATIONS[ll,jj,ii] <- Ybar[ii]+crossprod(Bb,CXXNEW[ii,])
            } 
          }
          SUMOFSQUARES[jj,ll] <- sum((Response-ESTIMATIONS[ll,jj,])^2)
        }
      }
    }
    nr <- nrow(SUMOFSQUARES)
    Numcolofmax <- max.col(1/SUMOFSQUARES)
    selectedpoint <- which.min(SUMOFSQUARES[cbind(1:nr,Numcolofmax)])
    Selected.points <- c(Selected.points,selectedpoint)
    Gridpoints <- Gridnew[-selectedpoint]
    XX <- STDCURVES[,Gridpoints] 
    XX.ALT.ROW <- XX[Ind11,]
    XX.REP.ROW <- XX[Ind12,]
    D2 <- (XX.ALT.ROW-XX.REP.ROW)^2
    Prox <- sqrt(apply(as.matrix(D2),1,sum))
    PROX <- matrix(Prox, nrow=nind)
    Max <- apply(PROX,1,max) 
    diag(PROX) <- Max 
    PROXSORTED <- apply(PROX,1,sort)
    Bw <- apply(PROXSORTED[Nbknn,],1,max)
    Bandwidths[kk] <- Bw[selectedpoint]
    Leaveoneout.estimated.values[[kk]] <- ESTIMATIONS[Numcolofmax[selectedpoint],selectedpoint,]
    Bcv[kk] <- sum((Response-Leaveoneout.estimated.values[[kk]])^2)/nind
    Bpcv[kk+1] <- Bcv[kk]*(1+(lforwardselection-kk)*pcvpar/log(nind))
    if((kk<lforwardselection) && (Bpcv[kk]<Bpcv[kk+1])) {
      exitloop <- T
      cat(paste("No!"),sep="\n")
      break      
    } else {
      cat(paste("Yes: design point",Gridnew[selectedpoint], "deleted"),sep="\n")
    }
  }
  if(exitloop) {
    Backwardselection <- Gridnew
    if(kk==1) {
      Bbwstd <- Fbwstd
      Bleaveoneout.estimated.values <- Fleaveoneout.estimated.values
    } else {
      Bbwstd <- Bandwidths[kk-1]
      Bleaveoneout.estimated.values <- Leaveoneout.estimated.values[[kk-1]]
    }
  } else {
    Backwardselection <- Gridpoints
    Bbwstd <- Bandwidths[kk]
    Bleaveoneout.estimated.values <- Leaveoneout.estimated.values[[kk]]
  }
  cat("Selected most-predicitve design points:", sep="\n")
  cat(Backwardselection, sep="\n")
  object <- list()
  class(object) <- "mpdp"
  object$Bwdselection <- Backwardselection
  object$Fwdselection <- Forwardselection
  object$bandwidth <- Bbwstd
  object$Bwd.estimated.values <- Bleaveoneout.estimated.values
  object$Bwdpcv <- Bpcv
  object$Fwdpcv <- Fpcv
  object$fpcv.status <- minimum.reached.for.pcv
  object$predictors <- CURVES
  object$responses <- Response
  object$kind.of.kernel <- kind.of.kernel
  return(object) 
}
#####################################################################
#####################################################################
#####################################################################
median.npfda <- function(SEMIMETRIC)
{
  Sum.of.distance <- apply(SEMIMETRIC, 1, sum)
  return(order(Sum.of.distance)[1])
}
#####################################################################
#####################################################################
#####################################################################
mplsr <- function(X, Y, K = 5)
{
  # Copyright (c) October 1993, Mike Denham.
  # Comments and Complaints to: snsdenhm@reading.ac.uk
  #
  # Orthogonal Scores Algorithm for PLS (Martens and Naes, pp. 121--123)
  #
  # X: predictors (matrix) 
  #
  # Y: multivariate response (matrix)
  #
  # K: The number of PLS factors in the model which must be less than or
  #    equal to the  rank of X.
  #
  # Returned Value is the vector of PLS regression coefficients
  #
  tol <- 1e-10
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  dx <- dim(X)
  nbclass <- ncol(Y)
  xbar <- apply(X, 2, sum)/dx[1]
  ybar <- apply(Y, 2, sum)/dx[1]
  X0 <- X - outer(rep(1, dx[1]), xbar)
  Y0 <- Y - outer(rep(1, dx[1]), ybar)
  W <- matrix(0, dx[2], K)
  P <- matrix(0, dx[2], K)
  Q <- matrix(0, nbclass, K)
  sumofsquaresY <- apply(Y0^2, 2, sum)
  u <- Y0[, order(sumofsquaresY)[nbclass]]
  tee <- 0
  for(i in 1:K) {
    test <- 1 + tol
    while(test > tol) {
      w <- crossprod(X0, u)
      w <- w/sqrt(crossprod(w)[1])
      W[, i] <- w
      teenew <- X0 %*% w
      test <- sum((tee - teenew)^2)
      tee <- teenew
      cee <- crossprod(tee)[1]
      p <- crossprod(X0, (tee/cee))
      P[, i] <- p
      q <- crossprod(Y0, tee)[, 1]/cee
      u <- Y0 %*% q
      u <- u/crossprod(q)[1]
    }
    Q[, i] <- q
    X0 <- X0 - tee %*% t(p)
    Y0 <- Y0 - tee %*% t(q)
  }
  COEF <- W %*% solve(crossprod(P, W)) %*% t(Q)
  b0 <- ybar - t(COEF) %*% xbar
  list(b0 = b0, COEF = COEF)
}
######################################################################################
######################################################################################
######################################################################################
######################################################################################
nested.funopadi = function(Classes, CONTAMINATED_CURVES, semimetric = "deriv", nbbwsmooth = 10, bw.adjust = F, kfold = 5, ...)
{
  ################################################################
  # Performs nonparametric functional discrimination of a sample of contaminated curves when 
  # a categorical response is observed (supervised classification). 
  # Two nested kernel estimators are implemented : one for the denoising step
  # and one for the nonparametric functional discrimination. 
  # The estimating algorithm selects automatically two smoothing parameters 
  # (one for the denoising step and one for the nonparametric functional 
  # discrimination step) by minimizing the k-fold cross-validation (the 
  # original dataset is split into "kfold" blocks or subsamples).
  # Remark: the measurements where the contaminated curves are observed must be   
  #         systematically the same (and equispaced) for the whole sample
  #    "Classes" vector containing the categorical responses
  #              giving the group number for each curve in 
  #              the matrix CONTAMINATED_CURVES (if nbclass is the number of 
  #              groups, "Classes" contains numbers 1,2,...,nbclass)
  #    "CONTAMINATED_CURVES" matrix containing the contaminated curves 
  #                          (stored row by row) used for the learning step
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  #    "nbbwsmooth" size of the bandwidths set used for selecting the optimal
  #                 for the denoising step (= 20 by default)
  #    "bw.adjust" logical value; if bw.adjust = T then the bandwidths used 
  #                in the denoising step are proportionnal to the standard deviation
  #                derived from the predictors observed at the measurement 
  #                where the kernel smoother is computed (by default bw.adjust = F)
  #    "kfold" integer giving the number of folds in the k-fold cross validation
  # Returns a list containing:
  #    "method" character string giving the mehod type (i.e. "discrimination")
  #    "args" list containing called arguments necessary to compute prediction with the "predict.npfda" function
  #    "Estimated.classnumber" vector containing estimated class membership 
  #                            for each curve of "CONTAMINATED_CURVES"
  #    "DENOISED_CURVES" matrix containing the denoised curves (stored row by row)
  #    "knn" optimal number of nearest neighbours in the regression step
  #    "Bandwidths" vector translating the optimal number of nearest neighbours in terms of 
  #                 bandwidths for each curve in the matrix "CONTAMINATED_CURVES" (regression step)
  #    "bws" optimal bandwidth used in the denoisiing step
  #    "misclas" k-fold misclassification rate
  ################################################################
  args = list(Responses = Classes, CONTAMINATED_CURVES = CONTAMINATED_CURVES, semimetric = semimetric, bw.adjust = bw.adjust, param = list(...))
  Classes = as.vector(Classes)
  nlearn = length(Classes)
  ### split the sample into kfold blocks ###
  nblock = nlearn %/% kfold
  listofblocks = list()
  oldsample = 1:nlearn
  for(kk in 1:(kfold - 1)){
    listofblocks[[kk]] = sample(oldsample, nblock)
    newsample = setdiff(oldsample, listofblocks[[kk]])
    oldsample = newsample
  }
  Alreadyused = c(listofblocks, recursive = T)
  listofblocks[[kfold]] = (1:nlearn)[- Alreadyused]
  # set the vector Knearest contains the sequence of k-nearest neighbours 
  # used for computing the optimal bandwidth
  nstart = nlearn - length(listofblocks[[1]])
  step = ceiling(nstart / 100)
  if(step == 0) step = 1
  Knearest = seq(from = trunc(0.05 * nstart), to = trunc(0.25 * nstart), by = step)
  nbknearest = length(Knearest)
  kmax = max(Knearest)
  # Grid of measurements where the contaminated curves are observed
  # Equispaced in [0, 1] and balanced (i.e. the same grid for all contaminated curves)
  Grid = seq(0, 1, length = ncol(CONTAMINATED_CURVES))
  ngrid = length(Grid)
  #################
  # denoising step	
  #################
  ### set the sequence of bandwidths for denoising stage
  bwsmooth_min = 5 * max(diff(Grid)) # smaller bandwidth = 5 times the difference between two consecutive grid points
  bwsmooth_max = 0.5 * diff(range(Grid)) # larger bandwidth = 25% of the grid range
  Bwsmooth = seq(from = bwsmooth_min, to = bwsmooth_max, length = nbbwsmooth)		
  if(bw.adjust){
    StandardDeviation = sqrt(apply(CONTAMINATED_CURVES, 2, var))
  }else{
    StandardDeviation = rep(1, ncol(CONTAMINATED_CURVES))		
  }
  Bws = Bwsmooth / mean(StandardDeviation)
  nbbws = length(Bws)
  # denoise the contaminated curves with the sequence of bandwidths contained in "Bws"
  SMOOTHED_CURVES = matrix(0, nrow = nlearn * ngrid, ncol = nbbws)
  for(bb in 1:nbbws) SMOOTHED_CURVES[, bb] = as.vector(kernel_smoother(Bws[bb], StandardDeviation, CONTAMINATED_CURVES, Grid))
  nbfold = length(listofblocks)
  nbclass = max(Classes)
  ESTIMATED.LABELS = matrix(0, nrow = nlearn, ncol = nbknearest)
  BINARY = matrix(0, nlearn, nbclass)
  Indices.list = c(listofblocks, recursive = T)
  for(g in 1:nbclass) BINARY[, g] = as.numeric(Classes == g)
  KFOLDCV = matrix(0, nbbws, nbknearest)
  sm = get(paste("semimetric.", semimetric, sep = ""))
  ###############################################################
  # combine (plug-in) the denoising step with the regression step
  ###############################################################
  #KFOLDCV = foreach(bb = 1:nbbws, .combine = 'cbind', .export = "sm") %dopar% {
  KFOLDCV = matrix(0, nbknearest, nbbws)
  for(bb in 1:nbbws){
    CURVES = matrix(SMOOTHED_CURVES[, bb], nrow = nlearn)
    # computing proximities between curves
    if(semimetric == "mplsr"){
      SEMIMETRIC1 = sm(Classes, CURVES, CURVES, ...)
    }else{
      if(semimetric == "L2"){
        SEMIMETRIC1 = sm(CURVES, CURVES)
      }else{
        SEMIMETRIC1 = sm(CURVES, CURVES, ...)
      }
    }
    # computing functional nonparametric supervised classification
    count = 0
    for(fold in 1:nbfold){
      Outsample = listofblocks[[fold]]
      Insample = (1:nlearn)[- Outsample]
      SEMIMETRICIO = SEMIMETRIC1[Insample, Outsample]
      BINARY.in = BINARY[Insample, ]
      nout = length(Outsample)
      HAT.PROB = matrix(0, nrow = nbclass, ncol = nbknearest)
      for(i in 1:nout) {
        Norm.diff = SEMIMETRICIO[, i]
        # "norm.order" gives the sequence k_1, k_2,... such that
        # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
        Norm.order = order(Norm.diff)
        # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
        # dq(X_{j_{kamx+2}},X_i)
        zz = Norm.diff[Norm.order][2:(kmax + 2)]
        # Bandwidth[l-1] contains (dq(X_{j_l},X_i) + 
        # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
        z = zz[ - (kmax + 1)]
        Bandwidth = 0.5 * (zz[-1] + z)
        ZMAT = matrix(rep(z, kmax), nrow = kmax, byrow = T)
        UMAT = ZMAT/Bandwidth
        KMAT = 1 - UMAT^2
        KMAT[UMAT > 1] = 0
        KMAT[col(KMAT) > row(KMAT)] = 0
        Ind.curves = Norm.order[2:(kmax + 1)]
        for(g in 1:nbclass) {
          Ind.resp = BINARY.in[Ind.curves, g]
          YMAT = matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
          HAT.PROB[g,  ] = .rowSums(YMAT[Knearest,  ] * KMAT[Knearest,  ], nbknearest, kmax)
        }
        Kmatsumbyrow = .rowSums(KMAT[Knearest,  ], nbknearest, kmax)
        HAT.PROB = HAT.PROB / matrix(Kmatsumbyrow, nbclass, nbknearest, byrow = T)
        count =  count + 1
        ESTIMATED.LABELS[count, ] = max.col(t(HAT.PROB))
      }
    }
    # k-fold criterion
    KFOLDCV[, bb] = .colSums(ESTIMATED.LABELS != Classes[Indices.list], nlearn, nbknearest) / nlearn
  }
  ### computing optimal parameters and corresponding k-fold CV
  Position = which(KFOLDCV == min(KFOLDCV), arr.ind = TRUE)
  if(nrow(Position) == 1){
    pos_knearest = Position[1, 1]
    pos_bws = Position[1, 2]
  }else{
    pos_knearest = Position[2, 1]
    pos_bws = Position[2, 2]
  }
  knearest_opt = Knearest[pos_knearest]
  bws_opt = Bws[pos_bws]	
  ###  Estimating labels with data driven parameters
  # denoising contaminated curves with optimal smoothing parameter
  CURVES = kernel_smoother(bws_opt, StandardDeviation, CONTAMINATED_CURVES, Grid)
  # computing proximities between denoised contaminated curves
  if(semimetric == "mplsr"){
    SEMIMETRIC1 = sm(Classes, CURVES, CURVES, ...)
  }else{
    if(semimetric == "L2"){
      SEMIMETRIC1 = sm(CURVES, CURVES)
    }else{
      SEMIMETRIC1 = sm(CURVES, CURVES, ...)
    }
  }
  # translate k-nearest neighbours in terms of continuous bandwidth for each 
  # curve in the learning sample
  Bandwidth.opt = 0
  for(i in 1:nlearn) {
    Sm1i = SEMIMETRIC1[, i]
    Bandwidth.opt[i] = sum(sort(Sm1i)[knearest_opt:(knearest_opt + 1)]) * 0.5
  }
  # Estimation of labels
  ARGUMENT = t(t(SEMIMETRIC1) / Bandwidth.opt)
  KERNEL = 1 - ARGUMENT^2
  KERNEL[ARGUMENT > 1] = 0
  diag(KERNEL) = 0
  Denom = .colSums(KERNEL, nlearn, nlearn)
  ESTIMATED.PROB = matrix(0, nrow = nlearn, ncol = nbclass)
  for(g in 1:nbclass) {
    PROBKERNEL = KERNEL * BINARY[, g]
    ESTIMATED.PROB[, g] = .colSums(PROBKERNEL, nlearn, nlearn) / Denom
  }
  Estimated.labels = max.col(ESTIMATED.PROB)
  # k-fold misclassification rate
  misclas.kfoldcv = sum(Estimated.labels != Classes) / nlearn
  # outputs
  outputs = list(method = "discrimination", args = args, Estimated.classnumber = Estimated.labels, DENOISED_CURVES = CURVES,
                 knn = knearest_opt, Bandwidths = Bandwidth.opt, bws = Bwsmooth[pos_bws], misclas = misclas.kfoldcv)
  class(outputs) = "npfda"
  return(outputs)
}
######################################################################################
######################################################################################
######################################################################################
######################################################################################
nested.funopare <- function(Responses, CONTAMINATED_CURVES, semimetric = "deriv", nbbwsmooth = 10, bw.adjust = F, kfold = 5, ...)
{
  ################################################################
  # Performs nonparametric functional regression of a sample of contaminated curves when 
  # a categorical response is observed (supervised classification). 
  # Nested kernel estimators are implemented : one for the denoising step
  # and one for the nonparametric functional discrimination. 
  # The estimating algorithm selects automatically two smoothing parameters 
  # (one for the denoising step and one for the nonparametric functional 
  # discrimination step) by minimizing the k-fold cross-validation (the 
  # original dataset is split into "kfold" blocks or subsamples).
  # Remark: the measurements where the contaminated curves are observed must be   
  #         systematically the same (and equispaced) for the whole sample
  #    "Responses" vector containing the categorical responses
  #              giving the group number for each curve in 
  #              the matrix CONTAMINATED_CURVES (if nbclass is the number of 
  #              groups, "Responses" contains numbers 1,2,...,nbclass)
  #    "CONTAMINATED_CURVES" matrix containing the contaminated curves 
  #                          (stored row by row) used for the learning step
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  #    "nbbwsmooth" size of the bandwidths set used for selecting the optimal
  #                 for the denoising step (= 20 by default)
  #    "bw.adjust" logical value; if bw.adjust = T then the bandwidths used 
  #                in the denoising step are proportionnal to the standard deviation
  #                derived from the predictors observed at the measurement 
  #                where the kernel smoother is computed (by default bw.adjust = F)
  #    "kfold" integer giving the number of folds in the k-fold cross validation
  # Returns a list containing:
  #    "method" character string giving the mehod type (i.e. "regression")
  #    "args" list containing called arguments necessary to compute prediction with the "predict.npfda" function
  #    "Estimated.values" vector containing estimated responses 
  #                       for each curve of "CONTAMINATED_CURVES"
  #    "DENOISED_CURVES" matrix containing the denoised curves (stored row by row)
  #    "knn" optimal number of nearest neighbours in the regression step
  #    "Bandwidths" vector translating the optimal number of nearest neighbours in terms of 
  #                 bandwidths for each curve in the matrix "CONTAMINATED_CURVES" (regression step)
  #    "bws" optimal bandwidth used in the denoisiing step
  #    "rmse" k-fold relative mean squared errors
  ################################################################
  args = list(Responses = Responses, CONTAMINATED_CURVES = CONTAMINATED_CURVES, semimetric = semimetric, bw.adjust = bw.adjust, param = list(...))
  Responses = as.vector(Responses)
  nlearn = length(Responses)
  ### split the sample into kfold blocks ###
  nblock = nlearn %/% kfold
  listofblocks = list()
  oldsample = 1:nlearn
  for(kk in 1:(kfold - 1)){
    listofblocks[[kk]] = sample(oldsample, nblock)
    newsample = setdiff(oldsample, listofblocks[[kk]])
    oldsample = newsample
  }
  Alreadyused = c(listofblocks, recursive = T)
  listofblocks[[kfold]] = (1:nlearn)[- Alreadyused]
  # set the vector Knearest contains the sequence of k-nearest neighbours 
  # used for computing the optimal bandwidth
  nstart = nlearn - length(listofblocks[[1]])
  step = ceiling(nstart / 100)
  if(step == 0) step = 1
  Knearest = seq(from = trunc(0.05 * nstart), to = trunc(0.25 * nstart), by = step)
  nbknearest = length(Knearest)
  kmax = max(Knearest)
  # Grid of measurements where the contaminated curves are observed
  # Equispaced in [0, 1] and balanced (i.e. the same grid for all contaminated curves)
  Grid = seq(0, 1, length = ncol(CONTAMINATED_CURVES))
  ngrid = length(Grid)
  #################
  # denoising step	
  #################
  ### set the sequence of bandwidths for denoising stage
  bwsmooth_min = 5 * max(diff(Grid)) # smaller bandwidth = 5 times the difference between two consecutive grid points
  bwsmooth_max = 0.5 * diff(range(Grid)) # larger bandwidth = 25% of the grid range
  Bwsmooth = seq(from = bwsmooth_min, to = bwsmooth_max, length = nbbwsmooth)		
  if(bw.adjust){
    StandardDeviation = sqrt(apply(CONTAMINATED_CURVES, 2, var))
  }else{
    StandardDeviation = rep(1, ncol(CONTAMINATED_CURVES))		
  }
  Bws = Bwsmooth / mean(StandardDeviation)
  nbbws = length(Bws)
  # denoise the contaminated curves with the sequence of bandwidths contained in "Bws"
  SMOOTHED_CURVES = matrix(0, nrow = nlearn * ngrid, ncol = nbbws)
  for(bb in 1:nbbws) SMOOTHED_CURVES[, bb] = as.vector(kernel_smoother(Bws[bb], StandardDeviation, CONTAMINATED_CURVES, Grid))
  nbfold = length(listofblocks)
  ESTIMATED.RESPONSES = matrix(0, nrow = nlearn, ncol = nbknearest)
  Indices.list = c(listofblocks, recursive = T)
  KFOLDCV = matrix(0, nbbws, nbknearest)
  sm = get(paste("semimetric.", semimetric, sep = ""))
  ###############################################################
  # combine (plug-in) the denoising step with the regression step
  ###############################################################
  #KFOLDCV = foreach(bb = 1:nbbws, .combine = 'cbind', .export = "sm") %dopar% {
  KFOLDCV = matrix(0, nbknearest, nbbws)
  for(bb in 1:nbbws){
    CURVES = matrix(SMOOTHED_CURVES[, bb], nrow = nlearn)
    # computing proximities between curves
    if(semimetric == "L2"){
      SEMIMETRIC1 = sm(CURVES, CURVES)
    }else{
      SEMIMETRIC1 = sm(CURVES, CURVES, ...)
    }
    # computing functional nonparametric supervised classification
    count = 0
    for(fold in 1:nbfold){
      Outsample = listofblocks[[fold]]
      Insample = (1:nlearn)[- Outsample]
      Responses.in = Responses[Insample]
      SEMIMETRICIO = SEMIMETRIC1[Insample, Outsample]
      nout = length(Outsample)
      for(i in 1:nout) {
        Norm.diff = SEMIMETRICIO[, i]
        # "norm.order" gives the sequence k_1, k_2,... such that
        # d(X_{k_1},X_i) < d(X_{k_2},X_i) < ...
        Norm.order = order(Norm.diff)
        # "zz" contains d(X_{k_2},X_i), d(X_{k_3},X_i),..., 
        # d(X_{j_{kamx+2}},X_i)
        zz = Norm.diff[Norm.order][2:(kmax + 2)]
        # Bandwidth[l-1] contains (dq(X_{j_l},X_i) + 
        # d(X_{j_l},X_i))/2 for l=2,...,kmax+2
        z = zz[ - (kmax + 1)]
        Bandwidth = 0.5 * (zz[-1] + z)
        ZMAT = matrix(rep(z, kmax), nrow = kmax, byrow = T)
        UMAT = ZMAT/Bandwidth
        KMAT = 1 - UMAT^2
        KMAT[UMAT > 1] = 0
        KMAT[col(KMAT) > row(KMAT)] = 0
        Ind.curves = Norm.order[2:(kmax + 1)]
        Ind.resp = Responses.in[Ind.curves]
        YMAT = matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
        Num = .rowSums(YMAT[Knearest,  ] * KMAT[Knearest,  ], nbknearest, kmax)
        Denom = .rowSums(KMAT[Knearest,  ], nbknearest, kmax)
        count =  count + 1
        ESTIMATED.RESPONSES[count,  ] = Num / Denom
      }
    }
    # k-fold criterion
    KFOLDCV[, bb] = .colSums((ESTIMATED.RESPONSES - Responses[Indices.list])^2, nlearn, nbknearest)	
  }
  ### computing optimal parameters and corresponding k-fold CV
  Position = which(KFOLDCV == min(KFOLDCV), arr.ind = TRUE)
  if(nrow(Position) == 1){
    pos_knearest = Position[1, 1]
    pos_bws = Position[1, 2]
  }else{
    pos_knearest = Position[2, 1]
    pos_bws = Position[2, 2]
  }
  knearest_opt = Knearest[pos_knearest]
  bws_opt = Bws[pos_bws]	
  ###  Computing estimations with data driven parameters
  # denoising contaminated curves with optimal smoothing parameter
  CURVES = kernel_smoother(bws_opt, StandardDeviation, CONTAMINATED_CURVES, Grid)
  # computing proximities between denoised contaminated curves
  if(semimetric == "L2"){
    SEMIMETRIC1 = sm(CURVES, CURVES)
  }else{
    SEMIMETRIC1 = sm(CURVES, CURVES, ...)
  }
  # translate k-nearest neighbours in terms of continuous bandwidth for each 
  # curve in the learning sample
  Bandwidth.opt = 0
  for(i in 1:nlearn) {
    Sm1i = SEMIMETRIC1[, i]
    Bandwidth.opt[i] = sum(sort(Sm1i)[knearest_opt:(knearest_opt + 1)]) * 0.5
  }
  # Computing estimations
  ARGUMENT = t(t(SEMIMETRIC1) / Bandwidth.opt)
  KERNEL = 1 - ARGUMENT^2
  KERNEL[ARGUMENT > 1] = 0
  diag(KERNEL) = 0
  Denom = .colSums(KERNEL, nlearn, nlearn)
  Estimated.responses = .colSums(KERNEL * Responses, nlearn, nlearn) / Denom
  # k-fold relative mean squared errors
  rmse = sum((Estimated.responses - Responses)^2) / (nlearn * var(Responses))
  # outputs
  outputs = list(method = "regression", args = args, Estimated.values = Estimated.responses, DENOISED_CURVES = CURVES,
                 knn = knearest_opt, Bandwidths = Bandwidth.opt, bws = Bwsmooth[pos_bws], rmse = rmse)
  class(outputs) = "npfda"
  return(outputs)
}
######################################################################################
######################################################################################
######################################################################################
######################################################################################
nke <- function(Responses, CONTAMINATED_CURVES, semimetric = "deriv", nbbwsmooth = 10, bw.adjust = F, kfold = 5, method = "regression", ...)
{
  #########################################################################
  # Estimate nonparametrically the relationship between responses and 
  # contaminated curves. 
  # Nested kernel estimators are implemented : one for the denoising step
  # and one for estimating the relationship. 
  # The estimating algorithm selects automatically two smoothing parameters 
  # (one for the denoising step and one for estimating the relationship) by 
  # minimizing a k-fold cross-validation (the original dataset is split into
  # "kfold" blocks or subsamples).
  # Remark: the measurements where the contaminated curves are observed must 
  #         be systematically the same (and equispaced) for the whole sample
  #    "Responses" vector containing responses corresponding to contaminated 
  #                curves
  #    "CONTAMINATED_CURVES" matrix containing the contaminated curves 
  #                          (stored row by row) used for the learning step
  #    "..." arguments needed for the call of the function computing 
  #          the semi-metric between curves
  #    "semimetric" character string allowing to choose the function 
  #                 computing the semimetric;  you can select 
  #                 "deriv" (default), "fourier", "hshift", "mplsr", 
  #                 and "pca"
  #    "nbbwsmooth" size of the bandwidths set used for selecting the optimal
  #                 for the denoising step (= 20 by default)
  #    "bw.adjust" logical value; if bw.adjust = T then the bandwidths used 
  #                in the denoising step are proportionnal to the standard 
  #                deviation derived from the predictors observed at the 
  #                measurement where the kernel smoother is computed (by 
  #                default bw.adjust = F)
  #    "kfold" integer giving the number of folds in the k-fold cross validation
  #    "method" character string equals to "regression" (default) or 
  #             "discrimination"
  # Returns outputs from the function "nested.funopare" if method = "regression"
  # (i.e. regression setting); in the other case (i.e. method = "discrimination")
  # outputs from the function "nested.funopadi" are given.
  ################################################################
  if(method == "regression"){
    outputs = nested.funopare(Responses, CONTAMINATED_CURVES, semimetric, nbbwsmooth, bw.adjust, kfold, ...)
  }else{
    outputs = nested.funopadi(Responses, CONTAMINATED_CURVES, semimetric, nbbwsmooth, bw.adjust, kfold, ...)
  }
  return(outputs)
}
#####################################################################
#####################################################################
############### LOCAL LINEAR MULTIVARIATE PREDICTION ################
############### USING AN OBJECT OF CLASS "mpdp"      ################
#####################################################################
#####################################################################
predict.mpdp <- function(object, CURVPRED)
{
  if(is.vector(CURVPRED)) CURVPRED <- as.matrix(CURVPRED) 
  VECEST <- object$predictors[,object$Bwdselection]
  VECPRED <- CURVPRED[,object$Bwdselection]
  Meancurve <- apply(VECEST,2,mean)
  CCOVARIATES <- t(t(VECEST)-Meancurve)
  Stdev <- sqrt(apply(CCOVARIATES^2,2,mean))
  COVARIATE <- t(t(CCOVARIATES)/Stdev)
  CVECPRED <- t(t(VECPRED)-Meancurve)
  PRED <- t(t(CVECPRED)/Stdev)
  nind <- nrow(COVARIATE)
  kernel <- get(object$kind.of.kernel)
  nind2 <- nrow(PRED)
  Ind1 <- rep(1:nind2,nind)
  Ind2 <- rep(1:nind,rep(nind2,nind))
  COVAR.REP.ROW <- COVARIATE[Ind2,]
  PRED.ALT.ROW <- PRED[Ind1,]
  D2 <- (PRED.ALT.ROW-COVAR.REP.ROW)^2
  D2OVERBW <- D2/(object$bandwidth^2)
  Dist <- sqrt(apply(D2OVERBW,1,sum))
  Ker <- rep(0,nind2*nind)
  Ker[Dist>0] <- kernel(Dist[Dist>0])
  KERNEL <- matrix(Ker, nrow=nind2)
  Denom <- apply(KERNEL,1,sum)
  NUM <- KERNEL %*% COVARIATE
  XBAR <- NUM/Denom
  Ynum <- as.vector(KERNEL %*% object$responses) 
  Ybar <- Ynum/Denom
  TP1 <- crossprod(t(KERNEL)*object$responses,COVARIATE)
  TP2 <- XBAR*Ynum
  TP1M2 <- TP1-TP2
  Predictions <- 0
  for(ii in 1:nind2){
    CCOVAR <- t(t(COVARIATE)-XBAR[ii,])
    TP3 <- crossprod(CCOVAR,(CCOVAR*KERNEL[ii,]))
    if(abs(det(TP3))<1e-8){
      Predictions[ii] <- Ybar[ii]
    }else{
      Bb <- solve(TP3,TP1M2[ii,])
      Cpred <- t(t(PRED)-XBAR[ii,])
      Predictions[ii] <- Ybar[ii]+crossprod(Bb,Cpred[ii,]) 
    } 
  }
  return(Predictions) 
}    
######################################################################################
######################################################################################
######################################################################################
######################################################################################
predict.npfda <- function(object, NEW_CONTAMINATED_CURVES)
  ##########################################################################################
# Compute predictions for an object of class "npfda" when providing new contaminted curves
##########################################################################################
#    "NEW_CONTAMINATED_CURVES" matrix containing new contaminated curves (stored row by row) 
#                              used for predicting corresponding scalar responses
# Returns a vector containing the predicted values
##########################################################################################
{
  npred = nrow(NEW_CONTAMINATED_CURVES)
  p = ncol(NEW_CONTAMINATED_CURVES)
  sm = get(paste("semimetric.", object$args$semimetric, sep = ""))
  param = object$args$param
  nlearn = nrow(object$args$CONTAMINATED_CURVES)
  knn = object$knn
  # Grid of measurements where the new contaminated curves are observed
  # ==> equispaced in [0, 1] and balanced (i.e. the same grid for all contaminated curves)
  Grid = seq(0, 1, length = p)
  ngrid = length(Grid)
  # store optimal smoothing parameter for the denoising step
  if(object$args$bw.adjust){
    StandardDeviation = sqrt(apply(NEW_CONTAMINATED_CURVES, 2, var))
  }else{
    StandardDeviation = rep(1, p)		
  }
  bws = object$bws / mean(StandardDeviation)
  # denoising contaminated curves (in the learning sample) with optimal smoothing parameter
  SMOOTHED_CURVES = kernel_smoother(bws, StandardDeviation, object$args$CONTAMINATED_CURVES, Grid)
  # denoising new contaminated curves with optimal smoothing parameter
  SMOOTHED_NEW_CURVES = kernel_smoother(bws, StandardDeviation, NEW_CONTAMINATED_CURVES, Grid)
  # computing proximities between curves
  if(object$args$semimetric == "L2"){
    SEMIMETRIC = sm(SMOOTHED_CURVES, SMOOTHED_NEW_CURVES)
  }else{
    if(length(param) == 3){
      SEMIMETRIC = sm(SMOOTHED_CURVES, SMOOTHED_NEW_CURVES, param[[1]], param[[2]], param[[3]])
    }else{
      SEMIMETRIC = sm(SMOOTHED_CURVES, SMOOTHED_NEW_CURVES, param[[1]])	
    }
  }
  # translate k-nearest neighbours in terms of continuous bandwidth for each new contaminated curve
  Bw.knn = 0
  for(k in 1:npred) {
    Sm2k = SEMIMETRIC[, k]
    Bw.knn[k] = sum(sort(Sm2k)[knn:(knn + 1)]) * 0.5
  }
  # Compute predictions
  ARGUMENT = t(t(SEMIMETRIC) / Bw.knn)
  KERNEL = 1 - ARGUMENT^2
  KERNEL[ARGUMENT > 1] = 0
  Denom = .colSums(KERNEL, nlearn, npred)
  if(object$method == "discrimination"){
    # Predicted labels
    nbclass = max(object$args$Responses)
    BINARY = matrix(0, nlearn, nbclass)
    for(g in 1:nbclass) BINARY[, g] = as.numeric(object$args$Responses == g)
    PREDICTED.PROB = matrix(0, nrow = npred, ncol = nbclass)
    for(g in 1:nbclass) {
      PROBKERNEL = KERNEL * BINARY[, g]
      PREDICTED.PROB[, g] = .colSums(PROBKERNEL, nlearn, npred) / Denom
    }
    Predictions = max.col(PREDICTED.PROB)	
  }else{
    # Predicted responses
    Predictions = .colSums(KERNEL * res$args$Responses, nlearn, npred) / Denom		
  }
  # outputs
  return(Predictions)
}
#####################################################################
#####################################################################
#####################################################################
prob.curve <- function(i, SEMIMETRIC, Hseq)
{
  ############################################################
  # Returns the probability curve for 
  # the ith unit :
  # "i" is an integer between 1 and size of the sample curves
  # "SEMIMETRIC" is the semimetric matrix of the sample curves
  # "Hseq" contains the sequence of the bandwidths "h"
  ############################################################
  Prox <- sort(SEMIMETRIC[i,  - i])
  Sbph <- 0
  for(k in 1:length(Hseq))
    Sbph[k] <- length(Prox[Prox < Hseq[k]])
  return(Sbph/ncol(SEMIMETRIC))
}
#####################################################################
#####################################################################
#####################################################################
quadratic <- function(u)
{
  #  quadratic kernel
  1 - (u)^2
}
####################################
### ASYMMETRICAL QUADRATIC KERNEL  #
####################################
quadratic2 <- function(u) 
{
  u[u<0] <- 1
  u[u>1] <- 1 
  return(1 - (u)^2) 
}
#####################################################################
#####################################################################
#####################################################################
rank.minima <- function(f)
{
  ######################################
  # Returns the rank of the local 
  # minima of the numeric strictly 
  # positive discretized function "f":
  # "Indmin" contains the ranks
  # ("f[Indmin]" gives the local minima)
  ######################################
  df <- diff(f)
  df[df>=0] <- 0
  df[df<0] <- 1
  ddf <- diff(df)
  Indmin <- (1:length(df))[ddf==-1]+1
  return(Indmin)
}
#####################################################################
#####################################################################
#####################################################################
semimetric.deriv <- function(DATA1, DATA2, q, nknot, range.grid)
{
  ###############################################################
  # Computes a semimetric between curves based on their derivatives.
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "q" order of derivation
  #    "nknot" number of interior knots (needed for defining the B-spline basis)
  #    "range.grid" vector of length 2 containing the range of the grid at 
  #                 which the curve are evaluated (i.e. range of the 
  #                 discretization)
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  library(splines)
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  #####################################################################
  # B-spline approximation of the curves containing in DATASET :
  # -----------------------------------------------------------
  # "knot" and "x" allow to define the B-spline basis
  # "coef.mat1[, i]" corresponds to the B-spline expansion
  # of the discretized curve contained in DATASET[i, ]. 
  # The B-spline approximation of the curve contained in "DATA1[i, ]" 
  # is given by "Bspline %*% coef.mat1[, i]"
  #####################################################################
  p <- ncol(DATA1)
  a <- range.grid[1]
  b <- range.grid[2]
  x <- seq(a, b, length = p)
  order.Bspline <- q + 3
  nknotmax <- (p - order.Bspline - 1)%/%2
  if(nknot > nknotmax){
    stop(paste("give a number nknot smaller than ",nknotmax, " for avoiding ill-conditioned matrix"))
  }
  Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
  delta <- sort(c(rep(c(a, b), order.Bspline), Knot))
  Bspline <- splineDesign(delta, x, order.Bspline)
  Cmat <- crossprod(Bspline)
  Dmat1 <- crossprod(Bspline, t(DATA1))
  coef.mat1 <- symsolve(Cmat, Dmat1)
  #######################################################################
  # Numerical integration by the Gauss method :
  # -------------------------------------------
  # The objects ending by "gauss" allow us to compute numerically  
  # integrals by means the "Gauss method" (lx.gauss=6 ==> the computation 
  # of the integral is exact for polynom of degree less or equal to 11).
  #######################################################################
  point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 
                   0.2386191861, 0.6612093865, 0.9324695142)
  weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
  x.gauss <- 0.5 * ((b + a) + (b - a) * point.gauss)
  lx.gauss <- length(x.gauss)
  Bspline.deriv <- splineDesign(delta, x.gauss, order.Bspline, rep(q, lx.gauss))
  H <- t(Bspline.deriv) %*% (Bspline.deriv * (weight.gauss * 0.5 * (b - a)))
  eigH <- eigen(H, sym = T)
  eigH$values[eigH$values < 0] <- 0
  Hhalf <- t(eigH$vectors %*% (t(eigH$vectors) * sqrt(eigH$values)))
  COEF1 <- t(Hhalf %*% coef.mat1)
  if(twodatasets){
    Dmat2 <- crossprod(Bspline, t(DATA2))
    coef.mat2 <- symsolve(Cmat, Dmat2)
    COEF2 <- t(Hhalf %*% coef.mat2)
  } else {
    COEF2 <- COEF1
  }
  SEMIMETRIC <- 0
  nbasis <- nrow(H)
  for(f in 1:nbasis)
    SEMIMETRIC <- SEMIMETRIC + outer(COEF1[, f], COEF2[, f], "-")^2
  return(sqrt(SEMIMETRIC))
}
#####################################################################
#####################################################################
#####################################################################
semimetric.fourier <- function(DATA1, DATA2, nderiv, nbasis, Range.grid, period=NULL)
{
  ###############################################################
  # Computes a semimetric between curves based on their Fourier expansion.
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "nderiv" order of derivation
  #    "nbasis" size of the basis
  #    "Range.grid" vector of length 2 containing the range of the grid at 
  #                 which the curve are evaluated (i.e. range of the 
  #                 discretization)
  #    "period" allows to select the period for the fourier expansion
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  p <- ncol(DATA1)
  nbasismax <- (p - nbasis)%/%2
  if(nbasis > nbasismax){
    stop(paste("give a number nbasis smaller than ",nbasismax, " for avoiding ill-conditioned matrix"))
  }
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  a <- Range.grid[1]
  b <- Range.grid[2]
  Eval <- seq(a, b, length = p)
  #####################################################################
  # Fourier approximation of the curves containing in DATA1 :
  # -----------------------------------------------------------
  # "fourier" allows to define the Fourier basis
  # "COEF1[, i]" corresponds to the Fourier expansion
  # of the discretized curve contained in DATA1[i, ]. 
  # The Fourier approximation of the curve contained in "DATA1[i, ]" 
  # is given by "FOURIER %*% COEF1[, i]"
  #####################################################################
  if(is.null(period)) period <- b - a
  FOURIER <- fourier(Eval, nbasis, period)
  CMAT <- crossprod(FOURIER)
  DMAT1 <- crossprod(FOURIER, t(DATA1))
  COEF1 <- symsolve(CMAT, DMAT1)
  #######################################################################
  # Numerical integration by the Gauss method :
  # -------------------------------------------
  # The objects ending by "gauss" allow us to compute numerically  
  # integrals by means the "Gauss method" (Leval.gauss=6 ==> the computation 
  # of the integral is exact for polynom of degree less or equal to 11).
  #######################################################################
  Point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 
                   0.2386191861, 0.6612093865, 0.9324695142)
  Weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,
                    0.360761573, 0.1713244924)
  Eval.gauss <- 0.5 * (b - a) * (1 + Point.gauss)
  Leval.gauss <- length(Eval.gauss)
  FOURIER.DERIV <- fourier(Eval.gauss, nbasis, period, nderiv)
  H <- t(FOURIER.DERIV) %*% (FOURIER.DERIV * (Weight.gauss * 0.5 * (b - a
  )))
  eigH <- eigen(H, sym = T)
  eigH$values[eigH$values < 0] <- 0
  HALF <- t(eigH$vectors %*% (t(eigH$vectors) * sqrt(eigH$values)))
  COEF1 <- t(HALF %*% COEF1)
  if(twodatasets) {
    DMAT2 <- crossprod(FOURIER, t(DATA2))
    COEF2 <- symsolve(CMAT, DMAT2)
    COEF2 <- t(HALF %*% COEF2)
  }
  else {
    COEF2 <- COEF1
  }
  SEMIMETRIC <- 0
  for(f in 1:nbasis)
    SEMIMETRIC <- SEMIMETRIC + outer(COEF1[, f], COEF2[, f], "-")^2
  return(sqrt(SEMIMETRIC))
}
#####################################################################
#####################################################################
#####################################################################
semimetric.hshift <- function(DATA1, DATA2, grid)
{
  ###############################################################
  # Computes between curves a semimetric taking into account an 
  # horizontal shift effect.   
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "grid" vector which defines the grid (one can choose 1,2,...,nbgrid
  #           where nbgrid is the number of points of the discretization) 
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  n1 <- nrow(DATA1)
  if(twodatasets) n2 <- nrow(DATA2) else n2 <- n1
  SEMIMETRIC <- matrix(0, nrow=n1, ncol=n2)
  if(!twodatasets){
    for(i in 1:(n1-1)){
      for(j in (i+1):n2){
        SEMIMETRIC[i,j] <- hshift(DATA1[i,], DATA2[j,], grid)$dist
      }  
    }
    SEMIMETRIC <- SEMIMETRIC + t(SEMIMETRIC)
  }else{
    for(i in 1:n1){
      for(j in 1:n2){
        SEMIMETRIC[i,j] <- hshift(DATA1[i,], DATA2[j,], grid)$dist
      }  
    }
  }
  return(SEMIMETRIC)
}
#####################################################################
#####################################################################
#####################################################################
semimetric.L2 = function(DATA1, DATA2)
{
  ###############################################################
  # Computes  L2 distance between curves.
  ################################################################
  #    "DATA1" matrix contains a first set of curves stored row by row
  #    "DATA2" matrix contains a second set of curves stored row by row
  # Returns a "semimetric" matrix containing the semimetrics
  # between all pairs (curve1, curve2) with curve1 in DATA1 and curve2
  # in DATA2 
  ###############################################################
  SEMIMETRIC = matrix(0, nrow(DATA1), nrow(DATA2))
  for(ii in 1:nrow(DATA2)){
    DIFF = t(DATA1) - DATA2[ii, ]
    SEMIMETRIC[, ii] = apply(DIFF^2, 2, sum) / (ncol(DATA2) - 1)
  }
  return(sqrt(SEMIMETRIC))
}
#####################################################################
#####################################################################
#####################################################################
semimetric.mplsr <- function(Classes1, DATA1, DATA2, q)
{
  ###############################################################
  # Computes between curves a semimetric based on the partial least 
  # squares method.
  #    "Classe1" vector containing a categorical response which 
  #              corresponds to class number for units stored in DATA1
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "q" the retained number of factors
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  qmax <- ncol(DATA1)
  if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
  n1 <- nrow(DATA1)
  nbclass <- max(Classes1)
  BINARY1 <- matrix(0, nrow = n1, ncol = nbclass)
  for(g in 1:nbclass) {
    BINARY1[, g] <- as.numeric(Classes1 == g)
  }
  mplsr.res <- mplsr(DATA1, BINARY1, q)
  COMPONENT1 <- DATA1 %*% mplsr.res$COEF
  COMPONENT1 <- outer(rep(1, n1), as.vector(mplsr.res$b0)) + COMPONENT1
  if(twodatasets) {
    n2 <- nrow(DATA2)
    COMPONENT2 <- DATA2 %*% mplsr.res$COEF
    COMPONENT2 <- outer(rep(1, n2), as.vector(mplsr.res$b0)) + 
      COMPONENT2
  }
  else {
    COMPONENT2 <- COMPONENT1
  }
  SEMIMETRIC <- 0
  for(g in 1:nbclass)
    SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, g], COMPONENT2[, 
                                                                 g], "-")^2
  return(sqrt(SEMIMETRIC))
}
#####################################################################
#####################################################################
#####################################################################
semimetric.pca <- function(DATA1, DATA2, q)
{
  ###############################################################
  # Computes between curves a pca-type semimetric based on the
  # functional principal components analysis method.
  #    "DATA1" matrix containing a first set of curves stored row by row
  #    "DATA2" matrix containing a second set of curves stored row by row
  #    "q" the retained number of principal components
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  ###############################################################
  if(is.vector(DATA1)) DATA1 <- as.matrix(t(DATA1))
  if(is.vector(DATA2)) DATA2 <- as.matrix(t(DATA2))
  testfordim <- sum(dim(DATA1)==dim(DATA2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(DATA1==DATA2)!=prod(dim(DATA1))
  qmax <- ncol(DATA1)
  if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
  n <- nrow(DATA1)
  COVARIANCE <- t(DATA1) %*% DATA1/n
  EIGENVECTORS <- eigen(COVARIANCE, sym = T)$vectors[, 1:q]
  COMPONENT1 <- DATA1 %*% EIGENVECTORS
  if(twodatasets) {
    COMPONENT2 <- DATA2 %*% EIGENVECTORS
  }
  else {
    COMPONENT2 <- COMPONENT1
  }
  SEMIMETRIC <- 0
  for(qq in 1:q)
    SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, qq], COMPONENT2[, 
                                                                  qq], "-")^2
  return(sqrt(SEMIMETRIC))
}
#####################################################################
#####################################################################
#####################################################################
summary.npfda <- function(object)
  ################################################################
# Display a summary of an object of class "npfda"
################################################################
{
  if(object$method == "discrimination"){
    cat(paste("k-fold misclassification rate:", object$misclas, "\n\n"))
  }else{
    cat(paste("k-fold relative mean squared erros:", object$rmse, "\n\n"))
  }
  cat("Denoising step\n")
  cat(paste("==> smoothing parameters:", format(object$bws, scientific = TRUE), "\n\n"))
  if(object$method == "discrimination"){
    cat("Functional nonparametric discrimination step\n")
  }else{
    cat("Functional nonparametric regression step\n")
  }
  cat(paste("==> number of nearest neighbours", object$knn, "\n\n"))
}
#####################################################################
#####################################################################
#####################################################################
symsolve <- function(Asym, Bmat)
{
  #  Performed by J.O. Ramsay and available on its 
  #  website http://ego.psych.mcgill.ca/misc/fda which contains a lot 
  #  of functions for analyzing functional data in a different way:
  #   Solves the system ASYM X = BMAT for X where ASYM is symmetric
  #   Returns X   
  n <- ncol(Asym)
  if(max(abs(Asym - t(Asym)))/max(abs(Asym)) > 1e-10)
    stop("Argument not symmetric.")
  Lmat <- chol(Asym, T)
  if(attr(Lmat, "rank") < n)
    stop("Argument singular.")
  Lmatinv <- solve(Lmat[, order(attr(Lmat, "pivot"))])
  Xmat <- Lmatinv %*% t(Lmatinv) %*% Bmat
  return(Xmat)
}
#####################################################################
#####################################################################
#####################################################################
triangle <- function(u)
{
  #  triangle kernel
  1 - u
}
#####################################################################
#####################################################################
#####################################################################
unbal2equibal <- function(DATA, INSTANTS, Range.grid, lnewgrid)
{
  ###############################################################
  # Transforms an unbalanced functional dataset into an equally
  # spaced balanced via linear interpolation.
  #    "DATA" matrix containing a set of unbalanced curves 
  #            stored row by row
  #    "INSTANTS" matrix containing the set of design points 
  #               stored row by row; INSTANTS[i,j]=t_{i,j} 
  #               (J_{max}=max_i J_i). As soon as J_i<J_{max}, 
  #               the ith row of INSTANTS is completed with NA's.
  #    "Range.grid" vector of length 2 containing the range of 
  #                 the desired grid (Range.grid[1]=t_1 and 
  #                 Range.grid[2]=t_J).
  #    "lnewgrid" length of the desired equally spaced grid 
  #               t_1,t_2,...,t_J.
  # Returns a matrix containing the new balanced curves stored 
  # row by row. 
  ###############################################################
  xmin <- Range.grid[1]
  xmax <- Range.grid[2]
  n <- nrow(DATA)
  NEWDATA <- matrix(0, n, lnewgrid)
  Newgrid <- seq(xmin, xmax, length = lnewgrid)
  for(i in 1:n) {
    Instantsi <- INSTANTS[i,  ]
    Instantsi <- Instantsi[!is.na(Instantsi)]
    Datai <- DATA[i,  ]
    Datai <- Datai[!is.na(Datai)]
    NEWDATA[i,  ] <- approx(Instantsi, Datai, xout = Newgrid, rule
                            = 2, f = 0)$y
  }
  return(list(NEWDATA = NEWDATA, Newgrid = Newgrid))
}





