#' Post-hoc comparisons and simultaneous confidence intervals for 
#' contrasts in multivariate factorial repeated measure designs
#' 
#' 
#' @author Mona Modrzik
#' 
#' @param object A \code{RM} object.
#' @param alpha Specifies the significance level, default is 0.05.
#' @param iter The number of iterations used for Wild Bootstrap. Default is 10000.
#' @param contrast The contrast matrix of interest. It can either be "pairwise" or "user-defined". 
#' @param contmat If contrast = "user-defined", the contrast matrix must be specified here. Note that
#' its rows must sum to zero.
#' @param type If contrast is "pairwise", the type of the pairwise comparison must be specified here. 
#' Calculation is based on the contrMat function in package multcomp, see the corresponding help page 
#' for details on the types of contrasts available.
#' @param base An integer specifying which group is considered the baseline group 
#' for Dunnett contrasts, see \code{\link[multcomp]{contrMat}}.
#' @param factor_name Specifies the factor for which post-hoc analysis are requested. Use All for all
#' pairwise comparisons.
#' @param seed A random seed for the resampling procedure. Optional, for reproducibility of the bootstrap procedure.
#' 
#' @details The RM.MCTP() function computes p-values and confidence intervals for the chosen contrast with MCTP and with Wild Bootstrap.
#'    
#' @return P-values and simultaneous confidence intervals for the chosen contrasts. 
#'
#'@importFrom stats aggregate
#'
#'@export

RM.MCTP <- function(object, alpha = 0.05, iter = 10000, contrast, contmat = NULL, type = NULL, factor_name = NULL, base  = 1, seed = NULL){
  
  if (!is(object, "RM")){
    stop("RM.MCTP is only implemented for objects of class RM, for MANOVA use simCI.")
  }
  
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a number between 0 and 1.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!(contrast %in% c("pairwise", "user-defined"))){
    stop("contrast must be either pairwise or user-defined!")
  }
  
  within <- object$within
  formula <- object$input$formula
  data <- object$input$data
  subject <- object$input$subject
  no.subf <- object$plotting$no.subf
  fl <- object$plotting$fl
  
  prep.dat <- prepare.data(formula = formula, data = data,
                                       subject = subject, within = within)
  d <- prep.dat[["lev.sub"]]
  no.whole <- prep.dat[["no.whole"]]
  lev <- prep.dat[["levels"]]
  a <- 1
  if (no.whole == 0) {
    a <- 1
  } else {
    for (i in 1:no.whole){
      a <- a * length(lev[[i]])
    }
  }
  
  en <- prep.dat[["n"]]
  n_list <- as.numeric(en)
  N <- sum(n_list)
  hypo <- prep.dat[["hypo"]]
  fac <- prep.dat[["fac"]]
  names(hypo) <- fac
  name_vars <- prep.dat[["EF"]]
  betw <- name_vars[1:no.whole]
  means_allg <- aggregate(formula, data = data, FUN = mean)
  response_pos <- attr(terms(formula), "response")
  response_name <- all.vars(formula)[response_pos]
  mean_order <- means_allg[do.call(order, means_allg[name_vars]),]
  mean_vec <- mean_order[[response_name]]
  
  if (contrast == "pairwise"){
    
    if(is.null(type)){
      stop("Please specify the type of pairwise comparison.")
    }
    
    if(is.null(factor_name)){
      stop("Please specify the factor.")
    }
    
    GHy <- hypo[[factor_name]]
    
    if (factor_name == "All"){
      aa <- prod(fl)
      labelfac <- do.call(paste, c(mean_order[name_vars], sep = "_"))
      labelfac <- unique(labelfac)
      ra <- rep(1,aa)
      names(ra) <- labelfac
      M <- multcomp::contrMat(ra, type = type, base =base)
      contmat <- M
    } else{
      factor_name_einzel <- unlist(strsplit(factor_name, ":"))
      if (!all(factor_name_einzel %in% names(fl)) ) {
        stop("Unknown factor in factor_name,", factor_name, "is not possible.")
      }
      labelfac <- do.call(paste, c(mean_order[factor_name_einzel], sep = "_"))
      labelfac <- unique(labelfac)
      k <- prod(fl[factor_name_einzel])
      rk <- rep(1,k)
      names(rk) <- labelfac
      if (length(factor_name_einzel) == 1){
        M <- multcomp::contrMat(rk, type = type, base =base)
        contmat <- M %*% GHy
      } else{
        M <- multcomp::contrMat(rep(1,k), type = type, base =base)
        contmat <- M %*% GHy
      }
    }
  }
  
  if (contrast == "user-defined"){
    if (is.null(contmat)){
      stop("Please specify a contrast matrix.")
    }
    if (is.null(dim(contmat))){
      contmat <- t(as.matrix(contmat))
    } else {
      if (sum(rowSums(contmat) > 1e-15) != 0){
        stop("The rows of the contrast matrix must sum to zero!")
      }
    }
    if (ncol(contmat) != prod(fl)){
      stop(paste("The contrast matrix must have", prod(fl), "columns."))
    }
  }
  
  create_block_diag <- function(matrices) {
    total_rows <- sum(sapply(matrices, nrow))
    total_cols <- sum(sapply(matrices, ncol))
    
    block_diag <- matrix(0, nrow = total_rows, ncol = total_cols)
    
    row_offset <- 0
    col_offset <- 0
    for (mat in matrices) {
      n_row <- nrow(mat)
      n_col <- ncol(mat)
      block_diag[(row_offset + 1):(row_offset + n_row), (col_offset + 1):(col_offset + n_col)] <- mat
      row_offset <- row_offset + n_row
      col_offset <- col_offset + n_col
    }
    return(block_diag)
  }
  
  Yw <- prep.dat[["data"]]
  V_hat <- array(0, dim = c(d, d, a))
  for (i in 1:a) {
    V_hat[,,i] <- cov(Yw[[i]])
  }
  
  Sigma_hat <- N*create_block_diag(lapply(1:a, function(i) V_hat[,,i] / n_list[i]))
  
  H <- contmat
  q <- nrow(H)
  sigma_hat_ell_ell <- numeric(q) 
  
  T_A <- numeric(q)
  for (j in 1:q) {
    sigma_hat_ell_ell[j] <- H[j, ] %*% Sigma_hat %*% (H[j, ])
    T_A[j] <- sqrt(N) * (H[j, ] %*% mean_vec) / sqrt(sigma_hat_ell_ell[j])
  }
  
  #Bootstrap
  
  T_A_star <- matrix(0, nrow = q, ncol = iter)  # Teststatistic for Bootstrap
  
  for (i in 1:iter) {
    data_star <- array(0, dim = c(N, d))
    W <- matrix(2*rbinom(N, 1, 1/2) - 1, ncol=1)  
    start_index <- 1
    
    kombi <- unique(data[betw])
    
    for (u in 1:a) {
      
      kombi_loop <- kombi[u,, drop = FALSE]
      grp <- rep(TRUE, nrow(mean_order))
      for(col in betw){
        grp <- grp & mean_order[[col]] == kombi_loop[[col]]
      }
      
      g_dat <- mean_order[grp,]
      
      group_mean <- g_dat[[response_name]]
      
      # One way layout
      if (no.whole == 0){
        group_mean <- mean_order[[response_name]]
      }
      
      response_matrix <- Yw[[u]]
      
      end_index <- start_index + n_list[u] - 1  
      data_star[start_index:end_index, ] <- W[start_index:end_index] * response_matrix - outer(W[start_index:end_index], group_mean)
      
      start_index <- end_index + 1 
    }
    
    group_indices <- split(1:(N), rep(1:a, times = n_list))  
    group_data_splits <- lapply(group_indices, function(idx) data_star[idx, ]) 
    
    V_hat_star <- simplify2array(lapply(group_data_splits, cov))  # Covariance
    means_star <- unlist(lapply(group_data_splits, colMeans))  # Mean
    names(means_star) <- NULL
    
    Sigma_hat_star <- N*create_block_diag(lapply(1:a, function(u) V_hat_star[,,u] / n_list[u]))
    
    # Teststatistic T^*
    for (j in 1:q) {
      T_A_star[j, i] <- sqrt(N) * (H[j, ] %*% means_star) / sqrt(H[j, ] %*% Sigma_hat_star %*% H[j, ])
    }
  }
  
  max_T_A <- max(abs(T_A))
  max_T_A_starABS <- apply(abs(T_A_star), 2, max)
  c_starABS <- quantile(max_T_A_starABS, probs = 1 - alpha)
  
  # Calculate Correlationmatrix R
  correlation_matrix <- outer(1:q, 1:q, Vectorize(function(j, m) {
    (H[j, ] %*% Sigma_hat %*% H[m, ]) /
      sqrt(sigma_hat_ell_ell[j] * sigma_hat_ell_ell[m])
  }))
  
  z_alpha_both <- mvtnorm::qmvnorm(1 - alpha, mean = rep(0, q), sigma = correlation_matrix, tail = "both")$quantile
  
  # p-values
  p_val_boot <- sum(max_T_A <= max_T_A_starABS) / iter
  p_val_individual_boot <- numeric(q)
  mittelpkt <- numeric(q)
  for (j in 1:q) {
    p_val_individual_boot[j] <- sum(abs(T_A[j]) <= abs(T_A_star[j,])) / iter
    mittelpkt[j] <- (H[j, ] %*% mean_vec)
  }
  p_val_individual_boot <- ifelse(p_val_individual_boot == 0,"<0.001", p_val_individual_boot)
  
  p_val_mctp <- 1 - mvtnorm::pmvnorm(lower = - max_T_A, upper = max_T_A, mean = rep(0,q),sigma = correlation_matrix)
  
  p_vals_indi_mctp  <- numeric(q)
  for (j in 1:q) {
    p_vals_indi_mctp[j] <- 1 - mvtnorm::pmvnorm(lower = -abs(T_A[j]), upper = abs(T_A[j]), mean = rep(0,q), sigma = correlation_matrix)
  }
  p_vals_indi_mctp <- ifelse(p_vals_indi_mctp < 1e-4,"<0.001", round(p_vals_indi_mctp, 4))
  
  # confidence intervals
  conf_intervals_boot <- matrix(0, nrow = q, ncol = 2)
  for (j in 1:q) {
    conf_intervals_boot[j, ] <- c((H[j, ] %*% mean_vec) - c_starABS * sqrt(sigma_hat_ell_ell[j] / N), 
                                  (H[j, ] %*% mean_vec) + c_starABS * sqrt(sigma_hat_ell_ell[j] / N))
  }
  
  conf_intervals_mctp <- matrix(0, nrow = q, ncol = 2)
  
  for (j in 1:q) {
    conf_intervals_mctp[j, ] <- c((H[j, ] %*% mean_vec) - z_alpha_both * sqrt(sigma_hat_ell_ell[j] / N), 
                                  (H[j, ] %*% mean_vec) + z_alpha_both * sqrt(sigma_hat_ell_ell[j] / N))
  }
  # Contrast used
  if (contrast == "user-defined"){
    out.contr <- "user-defined"
  }
  if (contrast == "pairwise"){
    out.contr <- type
  }
  
  if (contrast == "pairwise"){
    res_classic <- data.frame("contrast" = rownames(M), "Lower" =conf_intervals_mctp[,1], "Upper" =conf_intervals_mctp[,2], "p-value" = p_vals_indi_mctp)
    res_bootstrap <- data.frame("contrast" = rownames(M), "Lower" =conf_intervals_boot[,1], "Upper" =conf_intervals_boot[,2],"p-value" = p_val_individual_boot)
  }
  if (contrast == "user-defined"){
    res_classic <- data.frame("Lower" =conf_intervals_mctp[,1], "Upper" =conf_intervals_mctp[,2], "p-value" = p_vals_indi_mctp)
    res_bootstrap <- data.frame("Lower" =conf_intervals_boot[,1], "Upper" =conf_intervals_boot[,2],"p-value" = p_val_individual_boot)
  }
  
  
  cat("\n", "#------ Call -----#", 
      "\n", "\n", "-", "Contrast: ", out.contr, 
      "\n", "-", "Confidence level:", (1 - alpha) * 100, "%", 
      "\n") 
  
  cat("\n", "#------Post-hoc comparisons: MCTP -----#", 
      "\n", "\n") 
  print(res_classic)
  
  cat("\n", "#------Post-hoc comparisons: Wild-Bootstrap -----#", 
      "\n", "\n") 
  print(res_bootstrap)
}
