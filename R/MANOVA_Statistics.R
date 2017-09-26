## function for calculating test statistics for MANOVA
MANOVA.Stat<- function(data, n, hypo_matrix, iter, alpha, resampling, n.groups, p, CPU, seed, nf){
  
  N <- sum(n)
  H <- hypo_matrix
  x <- data
  #---------------- useful matrices ---------------------#
  A <-  t(rep(1 / n[1], n[1])) %x% diag(p)
  for (ii in 2:length(n)){
    B <- t(rep(1 / n[ii], n[ii])) %x% diag(p)
    A <- magic::adiag(A, B)
  }
  # -----------------------------------------------------#
  
  means <- A %*% x
  
  V <- list(NA)
  n.temp <- cumsum(c(0, n))
  for (i in 1:n.groups){
    y <- matrix(x[(n.temp[i]*p+1):(n.temp[i+1]*p)], ncol = p, byrow = TRUE)
    V[[i]] <- 1 / n[i] * cov(y)
  }
  
  sigma_hat <- V[[1]]
  for (i in 2:n.groups){
    sigma_hat <- magic::adiag(sigma_hat, V[[i]])
  }
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
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    
    # MATS
    DP <- diag(SnP)*diag(p*n.groups)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    return(list(WTPS, Q_N_P))
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
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- N * t(meansP) %*% TP %*% meansP
    # MATS
    DP <- diag(SnP)*diag(p*n.groups)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    return(list(WTPS, Q_N_P))
  }
  
  time1 <- Sys.time()  
  cl <- makeCluster(CPU)
  
  if(resampling == "paramBS"){
    if(seed != 0){
    parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    WTPS <- parSapply(cl, 1:iter, FUN = PBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)
    ecdf_MATS <- ecdf(unlist(WTPS[2, ]))
    p_valueMATS <- 1 - ecdf_MATS(Q_N)
  } else if(resampling == "WildBS"){
    if(seed != 0){
    parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    WTPS <- parSapply(cl, 1:iter, FUN = WBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)    
    ecdf_MATS <- ecdf(unlist(WTPS[2, ]))
    p_valueMATS <- 1-ecdf_MATS(Q_N)
  }
  
  time <- Sys.time() - time1
  parallel::stopCluster(cl)
  
  #------------------------ p-values -------------------------------#
  p_valueWTS <- 1 - pchisq(abs(WTS), df = df_WTS)

  #--------------------- Quantiles -------------------------------#
  quant_WTS <- quantile(ecdf_WTPS, alpha)
  quant_MATS <- quantile(ecdf_MATS, alpha)
  
  #-------------------- Output ----------------------------------#
  quantiles <- c(quant_WTS, quant_MATS)
  WTS_out <- c(WTS, df_WTS, p_valueWTS)
  MATS_out <- Q_N
  WTPS_out <- c(p_valueWTPS, p_valueMATS)
  result <- list(WTS = WTS_out, WTPS = WTPS_out, MATS = MATS_out,
                 Cov = Sn, Mean = means, time = time, quantiles = quantiles)
  return(result)
}


#-----------------------------------------------------------------------------------------------------------------------------
#######################
# for wide format data
#######################

MANOVA.Stat.wide <- function(Y, n, hypo_matrix, iter, alpha, resampling, CPU, seed){
  
  N <- sum(n)
  H <- hypo_matrix
  a <- length(Y)     # number of factor level combinations
  p <- ncol(Y[[1]]$response)   # number of dimensions
  
 means <- matrix(NA, ncol = p, nrow = a)
 
  for(i in 1:a){
    means[i, ] <- colMeans(Y[[i]]$response)
  }
 means <- c(t(means))
 
  V <- list(NA)
  for (i in 1:a){
    V[[i]] <- 1 / n[i] * cov(Y[[i]]$response)
  }
  
  sigma_hat <- V[[1]]
  for (i in 2:a){
    sigma_hat <- magic::adiag(sigma_hat, V[[i]])
  }
  Sn <- N * sigma_hat
  
  # WTS
  T <- t(H) %*% MASS::ginv(H %*% Sn %*% t(H)) %*% H
  WTS <- N * t(means) %*% T %*% means
  df_WTS <- Matrix::rankMatrix(H)[[1]]
  
  # MATS
  D <- diag(Sn)*diag(p*a)
  Q_N <- N* t(means)%*%t(H)%*%MASS::ginv(H%*%D%*%t(H))%*%H%*%means
  
  #--------------------------------- parametric bootstrap ---------------------------#
  PBS <- function(i, ...){
    # calculate mvrnorm for each group
    XP <- list()
    meansP <- list()
    for (i in 1:a){
      XP[[i]] <- MASS::mvrnorm(n[i], mu = rep(0, p), Sigma = n[i]*V[[i]])
      meansP[[i]] <- colMeans(XP[[i]])
    }
    meansP <- unlist(meansP)
    
    VP <- list()
    for(i in 1:a){
      VP[[i]] <- 1 / n[i] * cov(XP[[i]])
    }
    
    sigma_hatP <- VP[[1]]
    for (i in 2:a){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    
    # MATS
    DP <- diag(SnP)*diag(p*a)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    return(list(WTPS, Q_N_P))
  }
  
  #---------------------------------- Wild bootstrap ---------------------------------#
  WBS <- function(i, ...){
    
    VP <- list(NA)
    xperm <- list(NA)
    meansP <- list()
    for (i in 1:a){
      epsi <- 2*rbinom(n[i], 1, 1/2)-1
      xperm[[i]] <- epsi*(Y[[i]]$response - colMeans(Y[[i]]$response))
      VP[[i]] <- 1 / n[i] * cov(xperm[[i]])
      meansP[[i]] <- colMeans(xperm[[i]])
    }
    meansP <- unlist(meansP)
    
    sigma_hatP <- VP[[1]]
    for (i in 2:a){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- N * t(meansP) %*% TP %*% meansP
    # MATS
    DP <- diag(SnP)*diag(p*a)
    Q_N_P <- N* t(meansP)%*%t(H)%*%MASS::ginv(H%*%DP%*%t(H))%*%H%*%meansP
    return(list(WTPS, Q_N_P))
  }
  
  time1 <- Sys.time()  
  cl <- makeCluster(CPU)
  
  if(resampling == "paramBS"){
    if(seed != 0){
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    WTPS <- parSapply(cl, 1:iter, FUN = PBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)
    ecdf_MATS <- ecdf(unlist(WTPS[2, ]))
    p_valueMATS <- 1 - ecdf_MATS(Q_N)
  } else if(resampling == "WildBS"){
    if(seed != 0){
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }
    WTPS <- parSapply(cl, 1:iter, FUN = WBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)    
    ecdf_MATS <- ecdf(unlist(WTPS[2, ]))
    p_valueMATS <- 1-ecdf_MATS(Q_N)
  }
  
  time <- Sys.time() - time1
  parallel::stopCluster(cl)
  
  #------------------------ p-values -------------------------------#
  p_valueWTS <- 1 - pchisq(abs(WTS), df = df_WTS)
  #--------------------- Quantiles -------------------------------#
  quant_WTS <- quantile(ecdf_WTPS, alpha)
  quant_MATS <- quantile(ecdf_MATS, alpha)
  
  #-------------------- Output ----------------------------------#
  quantiles <- c(quant_WTS, quant_MATS)
  WTS_out <- c(WTS, df_WTS, p_valueWTS)
  MATS_out <- Q_N
  WTPS_out <- c(p_valueWTPS, p_valueMATS)
  result <- list(WTS = WTS_out, WTPS = WTPS_out, MATS = MATS_out,
                 Cov = Sn, Mean = means, time = time, quantiles = quantiles)
  return(result)
}