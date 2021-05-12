## function for calculating test statistics, permutation etc for repeated measures
RM.Stat<- function(Y, nind, hypo_matrix, iter, alpha, iii, hypo_counter, n.sub, 
                  resampling, para, CPU, seed, CI.method){
  
  N <- sum(nind)
  H <- hypo_matrix
  a <- length(Y)
  n <- rep(nind, each = n.sub)
  
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
  
  
  #---------------- useful matrices ---------------------#
   A <- diag(n.sub) %x% t(rep(1 / nind[1], nind[1]))
   if(a != 1){
   for (ii in 2:length(nind)){
     B <- diag(n.sub) %x% t(rep(1 / nind[ii], nind[ii]))
     A <- magic::adiag(A, B)
   } 
   }
  # -----------------------------------------------------#
  # WTS
  T <- t(H) %*% MASS::ginv(H %*% Sn %*% t(H)) %*% H
  WTS <- N * t(means) %*% T %*% means
  df_WTS <- Matrix::rankMatrix(H)[[1]]
  
  # ATS
  C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
  spur <- sum(diag(C %*% Sn))
  ATS <- N / spur * t(means) %*% C %*% means
  df_ATS <- spur ^ 2 / sum(diag(C %*% Sn %*% C %*% Sn))
   if (iii <= hypo_counter){
     # second df for independent factors
     Lambda <- diag(1 / (n - 1))
     D <- diag(C) * diag(ncol(C))
     df_ATS2 <- sum(diag(D %*% Sn)) ^ 2 / sum(diag(D %*% D %*% Sn %*% Sn %*% Lambda))
   } else {
    df_ATS2 <- Inf
  }
  
  #----------------------------Permutation --------------------------------#
  Perm <- function(arg, ...){
    x <- unlist(Y)
    xperm <- sample(x, replace = FALSE)
    meansP <- A %*% xperm
    n.temp <- n.sub*cumsum(c(0, nind))
    
    VP <- list()
    for(i in 1:a){
      yperm <- matrix(xperm[(n.temp[i]+1):n.temp[i+1]], ncol = n.sub)
      VP[[i]] <- cov(yperm)
    }
    
    sigma_hatP <- 1/nind[1]*VP[[1]]
    if(a != 1){
      for (i in 2:a){
        sigma_hatP <- magic::adiag(sigma_hatP, 1/nind[i]*VP[[i]])
      }
    }
    SnP <- N * sigma_hatP
    # WTPS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- N * t(meansP) %*% TP %*% meansP
    permout <- list(WTPS = WTPS, ATS_res = NA)
    return(permout)
  }
  
  #--------------------------------- parametric bootstrap ---------------------------#
  PBS <- function(ii, ...){
    # calculate mvrnorm for each group
    XP <- list()
    meansP <- list()
    for (i in 1:a){
      XP[[i]] <- MASS::mvrnorm(nind[i], mu = rep(0, n.sub), Sigma = V[[i]])
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
    
    # ATS
    C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
    spur <- sum(diag(C %*% SnP))
    ATS_res <- N / spur * t(meansP) %*% C %*% meansP
    
    pbs <- list(WTPS = WTPS, ATS_res = ATS_res, meansP = meansP)
    return(pbs)
  }
  
  #---------------------------------- Wild bootstrap ---------------------------------#
  WBS <- function(ii, ...){
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
    SnP <- N * sigma_hat
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- N * t(meansP) %*% TP %*% meansP
    # ATS
    C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
    spur <- sum(diag(C %*% SnP))
    ATS_res <- N / spur * t(meansP) %*% C %*% meansP
    wildbs <- list(WTPS = WTPS, ATS_res = ATS_res)
    return(wildbs)
  }
  
  if(para){
    cl <- makeCluster(CPU)
  
  if(seed != 0){
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  
  if(resampling == "Perm"){
    bs_out <- parallel::parLapply(cl, 1:iter, Perm)
  } else if(resampling == "paramBS"){
    bs_out <- parallel::parLapply(cl, 1:iter, PBS)
  } else if(resampling == "WildBS"){
    bs_out <- parallel::parLapply(cl, 1:iter, WBS)
  }
  WTPS <- parallel::parSapply(cl, bs_out, function(x) x$WTPS)
  ATSbs <- parallel::parSapply(cl, bs_out, function(x) x$ATS_res)
  parallel::stopCluster(cl)
  } else {
    set.seed(seed)
    if(resampling == "Perm"){
      bs_out <-lapply(1:iter, Perm)
    } else if(resampling == "paramBS"){
      bs_out <- lapply(1:iter, PBS)
    } else if(resampling == "WildBS"){
      bs_out <- lapply(1:iter, WBS)
    }
    WTPS <- sapply(bs_out, function(x) x$WTPS)
    ATSbs <- sapply(bs_out, function(x) x$ATS_res)
  }
  
  ecdf_WTPS <- ecdf(WTPS)
  p_valueWTPS <- 1-ecdf_WTPS(WTS)
  if(!is.na(ATSbs[1])){
    ecdf_ATS <- ecdf(ATSbs)
    p_valueATS_res <- 1 - ecdf_ATS(ATS)
  } else {
    p_valueATS_res <- NA
  }
  
  
  
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
