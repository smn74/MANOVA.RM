## function for calculating test statistics for MANOVA
MANOVA.Stat <- function(data, n, hypo_matrix, iter, alpha, resampling, n.groups, p,
                        para, CPU, seed){
  
  N <- sum(n)
  H <- hypo_matrix
  x <- data
  #---------------- useful matrices ---------------------#
  A <- Reduce(magic::adiag, lapply(seq_along(n), function(ii) t(rep(1/n[ii], n[ii])) %x% diag(p)))
  # -----------------------------------------------------#
  
  means <- A %*% x
  
  V <- list(NA)
  n.temp <- cumsum(c(0, n))
  for (i in 1:n.groups){
    y <- matrix(x[(n.temp[i]*p+1):(n.temp[i+1]*p)], ncol = p, byrow = TRUE)
    V[[i]] <- 1 / n[i] * cov(y)
  }
  
  sigma_hat <- Reduce(magic::adiag, V)
  Sn <- N * sigma_hat
  
  # WTS
  T <- t(H) %*% MASS::ginv(H %*% Sn %*% t(H)) %*% H
  WTS <- N * t(means) %*% T %*% means
  df_WTS <- Matrix::rankMatrix(H)[[1]]
  
  # MATS
  D <- diag(Sn)*diag(p*n.groups)
  Q_N <- N* t(means)%*%t(H)%*%MASS::ginv(H%*%D%*%t(H))%*%H%*%means
  
  #--------------------------------- parametric bootstrap ---------------------------#
  PBS <- function(i, ...){
    # calculate mvrnorm for each group
    XP <- list()
    meansP <- list()
    for (i in 1:n.groups){
      XP[[i]] <- MASS::mvrnorm(n[i], mu = rep(0, p), Sigma = n[i]*V[[i]])
      meansP[[i]] <- colMeans(XP[[i]])
    }
    meansP <- unlist(meansP)
    
    VP <- list()
    for(i in 1:n.groups){
      VP[[i]] <- 1 / n[i] * cov(XP[[i]])
    }
    
    sigma_hatP <- Reduce(magic::adiag, VP)
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    
    # MATS
    DP <- diag(SnP)*diag(p*n.groups)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    pbs <- list(WTPS = WTPS, Q_N_P = Q_N_P, meansP = meansP, DP = DP)
    return(pbs)
  }
  
  #---------------------------------- Wild bootstrap ---------------------------------#
  WBS <- function(i, ...){
    
    VP <- list(NA)
    xperm <- list(NA)
    meansP <- list()
    for (i in 1:n.groups){
      y <- matrix(x[(n.temp[i]*p+1):(n.temp[i+1]*p)], ncol = p, byrow = TRUE)
      epsi <- 2*rbinom(n[i], 1, 1/2)-1
      xperm[[i]] <- epsi*(y - colMeans(y))
      VP[[i]] <- 1 / n[i] * cov(xperm[[i]])
      meansP[[i]] <- colMeans(xperm[[i]])
    }
    meansP <- unlist(meansP)
    
    sigma_hatP <- Reduce(magic::adiag, VP)
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- N * t(meansP) %*% TP %*% meansP
    # MATS
    DP <- diag(SnP)*diag(p*n.groups)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    wbs <- list(WTPS = WTPS, Q_N_P = Q_N_P, meansP = meansP, DP = DP)
    return(wbs)
  }
  
  if(para){
    cl <- parallel::makeCluster(CPU)
    
    if(!is.null(seed)){
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    if(resampling == "paramBS"){
      bs_out <- parallel::parLapply(cl, 1:iter, PBS)
    } else if(resampling == "WildBS"){
      bs_out <- parallel::parLapply(cl, 1:iter, WBS)
    }
    
    WTPS <- sapply(bs_out, function(x) x$WTPS)
    MATSbs <- sapply(bs_out, function(x) x$Q_N_P)
    BSmeans <- lapply(bs_out, function(x) x$meansP)
    BSVar <- lapply(bs_out, function(x) x$DP)
    parallel::stopCluster(cl)
  } else {
    if(!is.null(seed)){
      set.seed(seed)
    }
    if(resampling == "paramBS"){
      bs_out <- lapply(1:iter, PBS)
    } else if(resampling == "WildBS"){
      bs_out <- lapply(1:iter, WBS)
    }
    
    WTPS <- sapply(bs_out, function(x) x$WTPS)
    MATSbs <- sapply(bs_out, function(x) x$Q_N_P)
    BSmeans <- lapply(bs_out, function(x) x$meansP)
    BSVar <- lapply(bs_out, function(x) x$DP)
  }
  ecdf_WTPS <- ecdf(WTPS)
  p_valueWTPS <- 1-ecdf_WTPS(WTS)
  ecdf_MATS <- ecdf(MATSbs)
  p_valueMATS <- 1 - ecdf_MATS(Q_N)
  
  
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
                 Cov = Sn, Mean = means, quantiles = quantiles, BSmeans = BSmeans,
                 BSVar = BSVar)
  return(result)
}


#-----------------------------------------------------------------------------------------------------------------------------
#######################
# for wide format data
#######################

# MANOVA.Stat.wide is mathematically identical to multRM.Statistic with t = 1
# (one observation per time point, i.e. no repeated measures within groups).
MANOVA.Stat.wide <- function(Y, n, hypo_matrix, iter, alpha, resampling,
                             para, CPU, seed, p){
  multRM.Statistic(Y, n, hypo_matrix, iter, alpha, resampling,
                   para, CPU, seed, p, t = 1)
}
