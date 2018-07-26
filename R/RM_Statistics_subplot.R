## function for calculating test statistics, permutation etc for repeated measures 
# with only sub-plot factors involved
RM.Stat.sub<- function(data, nind, n, hypo_matrix, iter, alpha, n.sub, n.groups,
                       resampling, seed, CI.method){
  
N <- nind
H <- hypo_matrix
x <- data
n.sub.total <- n.groups * n.sub

#---------------- useful matrices ---------------------#
A <- diag(n.sub) %x% t(rep(1 / nind, nind))
for (ii in 2:n.groups){
  B <- diag(n.sub) %x% t(rep(1 / nind, nind))
  A <- magic::adiag(A, B)
} 
# -----------------------------------------------------#
means <- A %*% x
y <- matrix(x, ncol = n.sub.total)
Sn <- cov(y)
  
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
df_ATS2 <- Inf


if(resampling == "Perm"){
  #----------------------------Permutation matrix--------------------------------#
  
  Perm <- matrix(0, nrow = N * n.sub.total, ncol = iter)
  for (pp in 1:iter){
    Perm[, pp] <- sample(1:(N * n.sub.total))
  }
  
  if(seed != 0){
    set.seed(seed)}
  #---------------------Wald-Type for permuted data ------------------------------#
  WTPS <- sapply(1:iter, function(arg){
    xperm <- x[Perm[, arg]]
    meansP <- A %*% xperm
    yperm <- matrix(xperm, ncol = n.sub.total)
    SnP <- cov(yperm)
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
  })
  ecdf_WTPS <- ecdf(WTPS)
  p_valueWTPS <- 1-ecdf_WTPS(WTS)
  p_valueATS_res <- NA
} else if(resampling == "paramBS"){
  if(seed != 0){
    set.seed(seed)}
  #--------------------------------- parametric bootstrap ---------------------------#
  WTPS <- sapply(1:iter, function(i, ...){
    xperm <- c(mvrnorm(N, rep(0, n.sub.total), Sn))
    meansP <- A %*% xperm
    yperm <- matrix(xperm, ncol = n.sub.total)
    SnP <- cov(yperm)
    
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    # ATS
    C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
    D <- diag(C) * diag(ncol(C))
    spur <- sum(diag(C %*% SnP))
    Lambda <- diag(1 / (n - 1))
    ATS_res <- N / spur * t(meansP) %*% C %*% meansP
    return(list(WTPS, ATS_res))
  })
  ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
  p_valueWTPS <- 1-ecdf_WTPS(WTS)    
  ecdf_ATS_res <- ecdf(unlist(WTPS[2, ]))
  p_valueATS_res <- 1-ecdf_ATS_res(ATS)
} else if(resampling == "WildBS"){
  if(seed != 0){
    set.seed(seed)}
  #---------------------------------- Wild bootstrap ---------------------------------#
  WTPS <- sapply(1:iter, function(i, ...){
    means2 <- rep(means, N)
    epsi <- 2*rbinom(N, 1, 1/2)-1
    xperm <- rep(epsi, n.sub.total)*(x-means2)
    meansP <- A %*% xperm
    yperm <- matrix(xperm, ncol = n.sub.total)
    SnP <- cov(yperm)
    
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
  })
  ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
  p_valueWTPS <- 1-ecdf_WTPS(WTS)    
  ecdf_ATS_res <- ecdf(unlist(WTPS[2, ]))
  p_valueATS_res <- 1-ecdf_ATS_res(ATS)  
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
               Cov = Sn, Mean = means, CI = CI)
return(result)
}