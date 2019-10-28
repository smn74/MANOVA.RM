## function for calculating test statistics, permutation etc for repeated measures
# Y1: data (only response vectors) split according to all groups (i.e. including time points)
# Y2: data split only into whole-plot factor groups
# n: observations in all the groups
# nind: number of individuals

multRM.Statistic <- function(Y, nind, hypo_matrix, iter, alpha, resampling, CPU, seed, p, t){
  
  N <- sum(nind)
  H <- hypo_matrix
  a <- length(Y)     # number of factor level combinations
  
  means <- sapply(Y, colMeans)
  means <- c(means)
  
  V <- lapply(Y, cov)
  sigma_hat <- 1/nind[1]*V[[1]]
  if(a != 1){
    for (i in 2:a){
      sigma_hat <- magic::adiag(sigma_hat, 1/nind[i]*V[[i]])
    }
  }
  Sn <- N * sigma_hat
  
  # WTS
  T <- t(H) %*% MASS::ginv(H %*% Sn %*% t(H)) %*% H
  WTS <- N * t(means) %*% T %*% means
  df_WTS <- Matrix::rankMatrix(H)[[1]]
  
  # MATS
  D <- diag(Sn)*diag(p*a*t)
  Q_N <- N* t(means)%*%t(H)%*%MASS::ginv(H%*%D%*%t(H))%*%H%*%means
  ### bootstrap functions
  
  #--------------------------------- parametric bootstrap ---------------------------#
  PBS <- function(i, ...){
    # calculate mvrnorm for each group
    XP <- list()
    meansP <- list()
    for (i in 1:a){
      XP[[i]] <- MASS::mvrnorm(nind[i], mu = rep(0, p*t), Sigma = V[[i]])
      meansP[[i]] <- colMeans(XP[[i]])
    }
    meansP <- unlist(meansP)
    
    VP <- lapply(XP, cov)
    
    sigma_hatP <- 1/nind[1]*VP[[1]]
    if (a != 1){
      for (i in 2:a){
        sigma_hatP <- magic::adiag(sigma_hatP, 1/nind[i]*VP[[i]])
      }
    }
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    
    # MATS
    DP <- diag(SnP)*diag(p*a*t)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    
    pbs <- list(WTPS = WTPS, Q_N_P = Q_N_P, meansP = meansP, DP = DP)
    return(pbs)
  }
  
  #---------------------------------- Wild bootstrap ---------------------------------#
  WBS <- function(i, ...){
    
    VP <- list(NA)
    xperm <- list(NA)
    meansP <- list()
    for (i in 1:a){
      epsi <- 2*rbinom(nind[i], 1, 1/2)-1
      xperm[[i]] <- epsi*(Y[[i]] - colMeans(Y[[i]]))
      VP[[i]] <- 1 / nind[i] * cov(xperm[[i]])
      meansP[[i]] <- colMeans(xperm[[i]])
    }
    meansP <- unlist(meansP)
    
    sigma_hatP <- VP[[1]]
    if (a != 1){
      for (i in 2:a){
        sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
      }
    }
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- N * t(meansP) %*% TP %*% meansP
    # MATS
    DP <- diag(SnP)*diag(p*a*t)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    wildbs <- list(WTPS=WTPS, Q_N_P=Q_N_P, meansP=meansP, DP=DP) 
    return(wildbs)
  }
  #-------------------------------------------------------------
  
  
  cl <- makeCluster(CPU)
  if(seed != 0){
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  if(resampling == "paramBS"){
    bs_out <- parallel::parLapply(cl, 1:iter, PBS)
  } else if(resampling == "WildBS"){
    bs_out <- parallel::parLapply(cl, 1:iter, WBS)
  }
 
  WTPS <- parallel::parSapply(cl, bs_out, function(x) x$WTPS)
  MATSbs <- parallel::parSapply(cl, bs_out, function(x) x$Q_N_P)
  ecdf_WTPS <- ecdf(WTPS)
  p_valueWTPS <- 1-ecdf_WTPS(WTS)
  ecdf_MATS <- ecdf(MATSbs)
  p_valueMATS <- 1 - ecdf_MATS(Q_N)
  
  BSmeans <- parallel::parLapply(cl, bs_out, function(x) x$meansP)
  BSVar <- parallel::parLapply(cl, bs_out, function(x) x$DP)
  parallel::stopCluster(cl)
  
  #------------------------ p-values -------------------------------#
  p_valueWTS <- 1 - pchisq(abs(WTS), df = df_WTS)
  #--------------------- Quantiles -------------------------------#
  quant_WTS <- quantile(ecdf_WTPS, 1-alpha)
  quant_MATS <- quantile(ecdf_MATS, 1-alpha)
  
  #-------------------- Output ----------------------------------#
  quantiles <- c(quant_WTS, quant_MATS)
  WTS_out <- c(WTS, df_WTS, p_valueWTS)
  MATS_out <- Q_N
  WTPS_out <- c(p_valueWTPS, p_valueMATS)
  result <- list(WTS = WTS_out, WTPS = WTPS_out, MATS = MATS_out,
                 Cov = Sn, Mean = means, time = time, quantiles = quantiles, BSmeans = BSmeans,
                 BSVar = BSVar)
  return(result)
}