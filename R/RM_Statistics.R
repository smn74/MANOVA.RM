## function for calculating test statistics, permutation etc for repeated measures
RM.Stat<- function(data, nind, n, hypo_matrix, iter, alpha, iii, hypo_counter, n.sub, 
                   n.groups, resampling, CPU, seed, CI.method){
  
  N <- sum(nind)
  H <- hypo_matrix
  x <- data
  
  #---------------- useful matrices ---------------------#
  A <- diag(n.sub) %x% t(rep(1 / nind[1], nind[1]))
  for (ii in 2:length(nind)){
    B <- diag(n.sub) %x% t(rep(1 / nind[ii], nind[ii]))
    A <- magic::adiag(A, B)
  } 
  # -----------------------------------------------------#
  means <- A %*% x
  
  V <- list(NA)
  n.temp <- n.sub*cumsum(c(0, nind))
  for (i in 1:n.groups){
    y <- matrix(x[(n.temp[i]+1):n.temp[i+1]], ncol = n.sub)
    V[[i]] <- 1 / nind[i] * cov(y)
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
  
  # ATS
  C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
  D <- diag(C) * diag(ncol(C))
  spur <- sum(diag(C %*% Sn))
  Lambda <- diag(1 / (n - 1))
  ATS <- N / spur * t(means) %*% C %*% means
  df_ATS <- spur ^ 2 / sum(diag(C %*% Sn %*% C %*% Sn))
  if (iii <= hypo_counter){
    df_ATS2 <- spur ^ 2 / sum(diag(D %*% D %*% Sn %*% Sn %*% Lambda))
  } else {
    df_ATS2 <- Inf
  }
  
  #----------------------------Permutation --------------------------------#
  Perm <- function(arg, ...){
    
    xperm <- sample(x, replace = FALSE)
    meansP <- A %*% xperm
    VP <- list(NA)
    for(i in 1:n.groups){
      yperm <- matrix(xperm[(n.temp[i]+1):n.temp[i+1]], ncol = n.sub)
      VP[[i]] <- 1 / nind[i] * cov(yperm)
    }
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    return(WTPS)
  }
  
  #--------------------------------- parametric bootstrap ---------------------------#
  PBS <- function(ii, ...){
    # calculate mvrnorm for each group
    XP <- list()
    meansP <- list()
    for (i in 1:n.groups){
      XP[[i]] <- MASS::mvrnorm(nind[i], mu = rep(0, n.sub), Sigma = nind[i]*V[[i]])
      meansP[[i]] <- colMeans(XP[[i]])
    }
    meansP <- unlist(meansP)
    
    VP <- list(NA)
    for(i in 1:n.groups){
      VP[[i]] <- 1 / nind[i] * cov(XP[[i]])
    }
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    # ATS
    C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
    D <- diag(C) * diag(ncol(C))
    spur <- sum(diag(C %*% SnP))
    Lambda <- diag(1 / (n - 1))
    ATS_res <- N / spur * t(meansP) %*% C %*% meansP
    return(list(WTPS, ATS_res))
  }
  
  #---------------------------------- Wild bootstrap ---------------------------------#
  WBS <- function(ii, ...){
    VP <- list(NA)
    xperm <- list(NA)
    for (i in 1:n.groups){
      y <- matrix(x[(n.temp[i]+1):n.temp[i+1]], ncol = n.sub)
      means2 <- rep(colMeans(y), nind[i])
      epsi <- 2*rbinom(nind[i], 1, 1/2)-1
      xperm[[i]] <- rep(epsi, n.sub)*(x[(n.temp[i]+1):n.temp[i+1]]-means2)
      yperm <- matrix(unlist(xperm[[i]])[(n.temp[i]+1):n.temp[i+1]], ncol = n.sub)
      VP[[i]] <- 1 / nind[i] * cov(yperm)
    }
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hat
    
    xperm <- unlist(xperm)
    meansP <- A %*% xperm
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    # ATS
    C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
    D <- diag(C) * diag(ncol(C))
    spur <- sum(diag(C %*% SnP))
    Lambda <- diag(1 / (n - 1))
    ATS_res <- N / spur * t(meansP) %*% C %*% meansP
    return(list(WTPS, ATS_res))
  }
  
  cl <- makeCluster(CPU)
  if(seed != 0){
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  
  if(resampling == "Perm"){
    WTPS <- parSapply(cl, 1:iter, FUN = Perm)
    ecdf_WTPS <- ecdf(WTPS)
    p_valueWTPS <- 1-ecdf_WTPS(WTS)
    p_valueATS_res <- NA
  } else if(resampling == "paramBS"){
    WTPS <- parSapply(cl, 1:iter, FUN = PBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)    
    ecdf_ATS_res <- ecdf(unlist(WTPS[2, ]))
    p_valueATS_res <- 1-ecdf_ATS_res(ATS)
  } else if(resampling == "WildBS"){
    WTPS <- parSapply(cl, 1:iter, FUN = WBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)    
    ecdf_ATS_res <- ecdf(unlist(WTPS[2, ]))
    p_valueATS_res <- 1-ecdf_ATS_res(ATS)    
  }
  
  parallel::stopCluster(cl)
  
  #------------------------ resampling quantile -------------------#
  quant_WTS <- quantile(ecdf_WTPS, 1-alpha)
  
  #------------------------ p-values -------------------------------#
  p_valueWTS <- 1 - pchisq(abs(WTS), df = df_WTS)
  p_valueATS <- 1 - pf(abs(ATS), df1 = df_ATS, df2 = df_ATS2)
  
  #---------------------- CIs -------------------------------------#
  if (CI.method == "t-quantile"){
    CI_lower <- means - sqrt(diag(Sn) / n) * qt(1 - alpha / 2, df = n)
    CI_upper <- means + sqrt(diag(Sn) / n) * qt(1 - alpha / 2, df = n)
  } else if (CI.method == "resampling"){
    CI_lower <- means - sqrt(diag(Sn) / n) * quant_WTS
    CI_upper <- means + sqrt(diag(Sn) / n) * quant_WTS
  }
  
  #-------------------- Output ----------------------------------#
  WTS_out <- c(WTS, df_WTS, p_valueWTS)
  ATS_out <- c(ATS, df_ATS, df_ATS2, p_valueATS)
  WTPS_out <- c(p_valueWTPS, p_valueATS_res)
  CI <- cbind(CI_lower, CI_upper)
  result <- list(WTS = WTS_out, WTPS = WTPS_out, ATS = ATS_out,
                 Cov = Sn, Mean = means, CI = CI, quantile = quant_WTS)
  return(result)
}