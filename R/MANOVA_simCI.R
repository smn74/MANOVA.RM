#' Multivariate post-hoc comparisons and simultaneous confidence intervals for contrasts in multivariate factorial designs
#' 
#' @param object A \code{MANOVA} object.
#' @param contrast The contrast matrix of interest, can either be "pairwise" or "user-defined". 
#' @param contmat If contrast = "user-defined", the contrast matrix must be specified here. Note that
#' its rows must sum to zero.
#' @param type If contrast is "pairwise", the type of the pairwise comparison must be specified here. 
#' Calculation is based on the contrMat function in package multcomp, see the corresponding help page 
#' for details on the types of contrasts available.
#' @param base An integer specifying which group is considered the baseline group 
#' for Dunnett contrasts, see \code{\link[multcomp]{contrMat}}.
#' @param interaction Logical. If interaction = FALSE in models with more than one factor, the factor of interest for 
#' the post-hoc analysis must be specified. Default is TRUE, which means post-hoc tests are performed for all factor
#' level combinations.
#' @param factor Only needed if interaction = FALSE. Specifies the factor for which post-hoc analysis are requested.
#' @param silent Set to TRUE to suppress output.
#' @param ... Not used yet.
#'   
#' @details The simCI() function computes the multivariate p-values for the chosen contrast of the multivariate mean vector 
#' based on the bootstrap version of the sum statistic. Details on this test can be found in Friedrich and Pauly (2018).
#' Furthermore, confidence intervals for summary effects (i.e., averaged over each dimension), also based on the bootstrap
#' version of the sum statistic, are returned as well.
#'    
#' @return Multivariate p-values and simultaneous confidence intervals for the chosen contrasts. 
#' 
#' @references Friedrich, S., and Pauly, M. (2018). MATS: Inference for potentially singular and
#'   heteroscedastic MANOVA. Journal of Multivariate Analysis, 165, 166-179.
#'   
#' @seealso \code{\link[multcomp]{contrMat}}
#' 
#' @importFrom multcomp contrMat
#' @importFrom methods is
#' @export
simCI <- function(object, contrast, contmat = NULL, type = NULL,
                  base = 1, interaction = TRUE, factor = NA, silent = FALSE, ...){
  
  if(!is(object, "MANOVA")){
    stop("simCI is currently only implemented for objects of class MANOVA.")
  }
  
  if(object$modeltype == "multRM"){
    stop("simCI is not yet implemented for repeated measures.")
  }
  
  if(object$nested){
    stop("The pairwise comparisons cannot be used in nested designs!")
  }
  
  if(!interaction){
    if(is.na(factor)){
      stop("Please enter the factor you wish to perform post-hoc tests on!")
    }
    if(!factor %in% object$factors){
      stop("The factor was not used in the original model!")
    }
    form1 <- strsplit(as.character(object$input$formula), "~", fixed = TRUE)[[2]]
    refit.formula <- as.formula(paste(form1, "~", factor))
    object <- object$modelcall(refit.formula, data = object$input$data, iter = object$input$iter, 
                               resampling = object$input$resampling, alpha = object$input$alpha,
                               seed = object$input$seed, subject = object$input$subject)
    
  }
  
  meanvec <- as.vector(t(object$Means))
  n <- object$Descriptive$n
  N <- sum(n)
  p <- object$p
  D <- diag(object$Covariance)*diag(p*length(n))
  factors <- object$factors
  nf <- object$nf
  BSmean <- object$BSMeans
  BSD <- object$BSVar
  alpha <- object$input$alpha
  fl <- object$fl
  lev <- subset(object$Descriptive, select = 1:nf)
  
  if(contrast == "user-defined"){
    if(is.null(contmat)){
      stop("Please specify a contrast matrix.")
    }
    # simple contrast vector
    if(is.null(dim(contmat))){
      contmat <- t(as.matrix(contmat))
    } else {
    if(sum(rowSums(contmat) > 1e-15) != 0){
      stop("The rows of the contrast matrix must sum to zero!")
    }
    }
    if (ncol(contmat) != p*prod(fl)){
      stop(paste("The contrast matrix must have", prod(fl)*p, "columns."))
    }
  }
  
  if(contrast == "pairwise"){
    if(is.null(type)){
      stop("Please specify the type of pairwise comparison, see the multcomp-package for
           details.")
    }
   # one-way
    if(nf == 1){
      names(n) <- lev[, 1]
    } else {
        names(n) <- do.call(paste, c(lev, sep = " "))
    }
      M <- contrMat(n, type = type, base)
      contmat <- M %x% t(rep(1, p))
  }
  
  # calculation of critical value based on sum statistic
  sumstat <- function(mean, var, ...){
    
    HDH <- diag(1/(contmat%*%var%*%t(contmat)))*diag(nrow(contmat))
    S_N_star <- N*t(contmat %*% mean)%*% HDH %*% contmat %*% mean 
    
    return(S_N_star)
  }
  
  S_star <- mapply(sumstat, BSmean, BSD)
  ecdf_S_star <- ecdf(S_star)
  q_star <- quantile(ecdf_S_star, 1-alpha)
  
  ## confidence intervals and multivariate p-values
  sci <- function(con, ...){
    center <- con%*%meanvec
    CI <- c(center - sqrt(q_star* t(con)%*%D%*%con/N), center + sqrt(q_star* t(con)%*%D%*%con/N))
    # test statistic for this contrast:
    Q_N_l <- N*t(center)%*%(t(con)%*%D%*%con)^(-1)%*%center
    p_val <- 1-ecdf_S_star(Q_N_l)
    
    result <- c(center, CI, p_val)
    names(result) <- c("Estimate", "Lower", "Upper", "p-value")
    return(result)
  }
  
  scis <- t(apply(contmat, 1, sci))
  if (contrast == "pairwise"){
    rownames(scis) <- rownames(M)
  }
  scis <- as.data.frame(scis)
  
  # output: split into multivariate p-values and confidence intervals separately
  p_val <- scis$`p-value`
  
  scis_out <- scis[, -4]
  
  contrast.output <- ifelse(contrast == "user-defined", "user-defined", type)
  if (contrast == "pairwise"){
    p_val<- data.frame("contrast" = rownames(M), "p-value" = p_val)
  }
  
  # avoid printing zeros
  # p_val[p_val[, "p.value"] == 0, "p.value"] <- "<0.001"
  
  if(!silent){
  cat("\n", "#------ Call -----#", 
      "\n", "\n", "-", "Contrast: ", contrast.output, 
      "\n", "-", "Confidence level:", (1 - alpha) * 100, "%", 
      "\n") 
  
  cat("\n", "#------Multivariate post-hoc comparisons: p-values -----#", 
  "\n", "\n") 
  print(p_val)
  
  cat("\n", "#-----------Confidence intervals for summary effects-------------#", 
  "\n", "\n")
  print(scis_out)
  }
  
 invisible(scis)
}
