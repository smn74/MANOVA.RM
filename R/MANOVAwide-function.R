#' Tests for Multivariate Data in Semi-Parametric Factorial Designs
#' 
#' The MANOVA.wide function calculates the Wald-type statistic (WTS)
#'  and a modified ANOVA-type statistic (MATS) as well as resampling versions of 
#' these test statistics for 
#' semi-parametric multivariate data provided in wide format.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side 
#'   contains the matrix of response variables and the right hand side contains the factor 
#'   variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. Data must be in wide format.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "paramBS"
#'   (parametric bootstrap approach) and "WildBS" (wild bootstrap approach with
#'   Rademacher weights). The Wild Bootstrap is calculated for all test statistics.
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
#' @param seed A random seed for the resampling procedure. If omitted, no 
#'   reproducible seed is set.
#' @param nested.levels.unique A logical specifying whether the levels of the nested factor(s)
#'   are labeled uniquely or not. Default is FALSE, i.e., the levels of the nested 
#'   factor are the same for each level of the main factor. For an example and more explanations
#'   see the GFD package and the corresponding vignette.
#' @param dec Number of decimals the results should be rounded to. Default is 3.
#'  
#' @section NOTE: The number of resampling iterations has been set to 100 in the examples due to run time 
#' restrictions on CRAN. Usually it is recommended to use at least 1000 iterations.
#'     
#' @return See \code{\link{MANOVA}}
#'  
#' @examples 
#' #Example on producing plastic film from Krzanowski (1998, p. 381), see \code{\link{manova.summary}}
#' tear <- c(6.5, 6.2, 5.8, 6.5, 6.5, 6.9, 7.2, 6.9, 6.1, 6.3,
#'           6.7, 6.6, 7.2, 7.1, 6.8, 7.1, 7.0, 7.2, 7.5, 7.6)
#' gloss <- c(9.5, 9.9, 9.6, 9.6, 9.2, 9.1, 10.0, 9.9, 9.5, 9.4,
#'            9.1, 9.3, 8.3, 8.4, 8.5, 9.2, 8.8, 9.7, 10.1, 9.2)
#' opacity <- c(4.4, 6.4, 3.0, 4.1, 0.8, 5.7, 2.0, 3.9, 1.9, 5.7,
#'              2.8, 4.1, 3.8, 1.6, 3.4, 8.4, 5.2, 6.9, 2.7, 1.9)
#' rate     <- gl(2,10, labels = c("Low", "High"))
#' additive <- gl(2, 5, length = 20, labels = c("Low", "High"))
#' example <- data.frame(tear, gloss, opacity, rate, additive)
#' fit <- MANOVA.wide(cbind(tear, gloss, opacity) ~ rate * additive, 
#' data = example, iter = 100, CPU = 1)
#' summary(fit)
#'
#' @seealso \code{\link{MANOVA}}
#'
#' @export

MANOVA.wide <- function(formula, data,
                   iter = 10000, alpha = 0.05, resampling = "paramBS", CPU,
                   seed, nested.levels.unique = FALSE, dec = 3){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  input_list <- list(formula = formula, data = data,
                     iter = iter, alpha = alpha, resampling = resampling)
  
  test1 <- hasArg(CPU)
  if(!test1){
    CPU <- parallel::detectCores()
  }
  
  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }
  
  dat <- model.frame(formula, data)
  nr_hypo <- attr(terms(formula), "factors")
  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  fac_names <- colnames(nr_hypo)
  
  outcome_names <- rownames(nr_hypo)[1]  # names of outcome variables
  # extract names of outcome variables
  split1 <- strsplit(outcome_names, "(", fixed = TRUE)[[1]][-1]
  split2 <- strsplit(split1, ")", fixed = TRUE)[[1]]
  split3 <- strsplit(split2, ",")[[1]]
  
  EF <- rownames(nr_hypo)[-1]  # names of influencing factors
  nf <- length(EF)
  names(dat) <- c("response", EF)
  #no. dimensions
  p <- ncol(dat$response)
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, (aa + 1)]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, (jj + 1)]))
  }
  lev_names <- expand.grid(levels)
  
  if (nf == 1) {
    # one-way layout
    dat2 <- dat[order(dat[, 2]), ]
    fac.groups <- dat2[, 2]
    n.groups <- prod(fl)
    Y <- split(dat2, fac.groups)
    n <- sapply(Y, nrow)
    hypo <- (diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(p)
    
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    MATS_out <- NA
    WTPS_out <- rep(NA, 2)
    quantiles <- matrix(NA, 2, 1)
    rownames(WTS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- MANOVA.Stat.wide(Y, n = n, hypo, iter = iter, alpha, resampling, CPU, seed)    
    WTS_out <- round(results$WTS, dec)
    MATS_out <- round(results$MATS, dec)
    WTPS_out <- round(results$WTPS, dec)
    quantiles <- results$quantiles
    names(quantiles) <- c("WTS_resampling", "MATS_resampling")
    mean_out <- matrix(round(results$Mean, dec), ncol = p, byrow = TRUE)
    Var_out <- results$Cov
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(EF, "n", split3)  
    names(WTS_out) <- cbind ("Test statistic", "df",
                             "p-value")
    names(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$Covariance <- Var_out
    output$WTS <- WTS_out
    output$MATS <- MATS_out
    output$resampling <- WTPS_out
    output$quantile <- quantiles
    output$nf <- nf
    output$factors <- fac_names
    output$H <- hypo
    output$p <- p
    output$fl <- fl
    output$Means <- mean_out
    # end one-way layout ------------------------------------------------------
  } else {
    dat2 <- dat[do.call(order, dat[, 2:(nf + 1)]), ]
    fac.groups <- do.call(list, dat2[, 2:(nf+1)])
   
    Y <- split(dat2, fac.groups, lex.order = TRUE)
    n <- sapply(Y, nrow)
    
    nested <- grepl(":", formula)
    
    if (sum(nested) > 0) {
      # nested
      
      # if nested factor is named uniquely
      if (nested.levels.unique){
        # delete factorcombinations which don't exist
        n <- n[n != 0]
        # create correct level combinations
        blev <- list()
        lev_names <- list()
        for (ii in 1:length(levels[[1]])) {
          blev[[ii]] <- levels(as.factor(dat[, 3][dat[, 2] == levels[[1]][ii]]))
          lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
        }
        if (nf == 2) {
          lev_names <- as.factor(unlist(lev_names))
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev)
        } else {
          lev_names <- lapply(lev_names, rep,
                              length(levels[[3]]) / length(levels[[2]]))
          lev_names <- lapply(lev_names, sort)
          lev_names <- as.factor(unlist(lev_names))
          blev <- lapply(blev, rep, length(levels[[3]]) / length(levels[[2]]))
          blev <- lapply(blev, sort)
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev, as.factor(levels[[3]]))
        }
        # correct for wrong counting of nested factors
        if (nf == 2) {
          fl[2] <- fl[2] / fl[1]
        } else if (nf == 3) {
          fl[3] <- fl[3] / fl[2]
          fl[2] <- fl[2] / fl[1]
        }
      }
      hypo_matrices <- HN_MANOVA(fl, p)
    } else {
      # crossed
      hypo_matrices <- HC_MANOVA(fl, perm_names, fac_names, p)[[1]]
    }
    
    n.groups <- prod(fl)
    
    if(length(Y) != n.groups){
      index <- NULL
      for(i in 1:length(Y)){
        if(nrow(Y[[i]]) == 0){
           index <- c(index, i)
        }
      }
      Y <- Y[-index]
    }
    
    
    # ---------------------- error detection ------------------------------------
    
    # mixture of nested and crossed designs is not possible
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is
           not implemented!")
    }
    # only 3-way nested designs are possible
    if (sum(nested) > 0 && nf >= 4) {
      stop("Four- and higher way nested designs are
           not implemented!")
    }
    # no factor combinations with less than 2 observations
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination
           with less than 2 observations!")
    }
    
    #--------------------------------------------------------------------------#
    
    WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
    WTPS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2)
    MATS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 1)
    quantiles <- matrix(NA, ncol = 2, nrow = length(hypo_matrices))
    rownames(WTS_out) <- fac_names
    rownames(WTPS_out) <- fac_names
    rownames(MATS_out) <- fac_names
    rownames(quantiles) <- fac_names
    colnames(MATS_out) <- "Test statistic"
    colnames(quantiles) <- c("WTS_resampling", "MATS_resampling")
    # calculate results
    for (i in 1:length(hypo_matrices)) {
      results <- MANOVA.Stat.wide(Y, n, hypo_matrices[[i]],
                             iter, alpha, resampling, CPU, seed)
      WTS_out[i, ] <- round(results$WTS, dec)
      WTPS_out[i, ] <- round(results$WTPS, dec)
      MATS_out[i] <- round(results$MATS, dec)
      quantiles[i, ] <- results$quantiles
    }
    # time needed for resampling calculations
    time <- results$time
    mean_out <- matrix(round(results$Mean, dec), ncol = p, byrow = TRUE)
    Var_out <- results$Cov
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(EF, "n", split3)
    rownames(descriptive) <- NULL
    
    # Output ------------------------------------------------------
    colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
    colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
    output <- list()
    output$time <- time
    output$input <- input_list
    output$Descriptive <- descriptive
    output$Covariance <- Var_out
    output$Means <- mean_out
    output$MATS <- MATS_out
    output$WTS <- WTS_out
    output$resampling <- WTPS_out
    output$quantile <- quantiles
    output$nf <- nf
    output$H <- hypo_matrices
    output$factors <- fac_names
    output$p <- p
  }
  class(output) <- "MANOVA"
  return(output)
}
