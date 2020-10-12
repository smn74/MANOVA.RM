#' Tests for Multivariate Data in Semi-Parametric Factorial Designs
#' 
#' The MANOVA function calculates the Wald-type statistic (WTS) and a modified ANOVA-type 
#' statistic (MATS) as well as resampling versions of 
#' these test statistics for 
#' semi-parametric multivariate data.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side 
#'   contains the response variable and the right hand side contains the factor 
#'   variables of interest.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. Data must be in long format and must not contain missing values.
#' @param subject The column name of the subjects in the data.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "paramBS"
#'   (parametric bootstrap approach) and "WildBS" (wild bootstrap approach with
#'   Rademacher weights).
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
#' @details The MANOVA() function provides the Wald-type statistic (WTS) as well as
#'   the modified ANOVA-type statistic (MATS) for multivariate designs with metric data as described in 
#'   Konietschke et al. (2015) and Friedrich and Pauly (2018), respectively. The MATS is invariant
#'   under scale transformations of the components and applicable to designs with singular covariance
#'   matrices.
#'   Both tests are applicable for non-normal error terms, 
#'   different sample sizes and/or heteroscedastic variances.
#'   They are implemented for designs with an arbitrary number of
#'   crossed factors or for nested designs. In addition to the asymptotic
#'   p-values, the function also provides p-values based on resampling approaches.
#'  
#' @section NOTE: The number of resampling iterations has been set to 10 in the examples due to run time 
#' restrictions on CRAN. Usually it is recommended to use at least 1000 iterations. 
#' For more information and detailed examples also refer to the package vignette.
#'     
#' @return A \code{MANOVA} object containing the following components: 
#'   \item{Descriptive}{Some descriptive statistics of the data for all factor 
#'   level combinations. Displayed are the number of individuals per factor 
#'   level combination and the vector of means (one column per dimension).}
#'   \item{Covariance}{The estimated covariance matrix.} 
#'   \item{WTS}{The value of the WTS along with degrees of freedom of the
#'   central chi-square distribution and p-value.} 
#'   \item{MATS}{The value of the MATS.} 
#'   \item{resampling}{p-values for the test statistic based on the
#'   chosen resampling approach.}
#'  
#' @examples data(EEG)
#' EEG_mod <- MANOVA(resp ~ sex * diagnosis, 
#'                     data = EEG, subject = "id", resampling = "paramBS", 
#'                     alpha = 0.05, iter = 10, CPU = 1)
#' summary(EEG_mod)
#' 
#' @seealso \code{\link{RM}}
#'   
#' @references
#' Friedrich, S., Konietschke, F., and Pauly, M. (2019). Resampling-Based Analysis of Multivariate Data 
#' and Repeated Measures Designs with the R Package MANOVA.RM. The R Journal, 11(2), 380-400.
#' 
#'   Konietschke, F., Bathke, A. C., Harrar, S. W. and Pauly, M. (2015). 
#'   Parametric and nonparametric bootstrap methods for general MANOVA. Journal 
#'   of Multivariate Analysis, 140, 291-301.
#'   
#'   Friedrich, S., Brunner, E. and Pauly, M. (2017). Permuting longitudinal data
#'   in spite of the dependencies. Journal of Multivariate Analysis, 153, 255-265.
#'   
#'   Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and 
#'   Hoeller, Y. (2018). Testing Mean Differences among Groups: Multivariate and Repeated 
#'   Measures Analysis with Minimal Assumptions. Multivariate Behavioral Research, 53(3), 348-359, 
#'   Doi: 10.1080/00273171.2018.1446320.
#'   
#'   Friedrich, S., Konietschke, F., Pauly, M. (2017). GFD - An 
#'   R-package for the Analysis of General Factorial Designs. 
#'   Journal of Statistical Software, 79(1), 1-18.
#'   
#'   Friedrich, S., and Pauly, M. (2018). MATS: Inference for potentially singular and
#'   heteroscedastic MANOVA. Journal of Multivariate Analysis, 165, 166-179.
#'    
#'   
#' @importFrom graphics axis legend par plot title abline points
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var cov rbinom quantile as.formula
#' @importFrom utils read.table
#' @importFrom methods hasArg
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster parSapply detectCores
#' @importFrom ellipse ellipse
#'   
#' @export

MANOVA <- function(formula, data, subject,
                   iter = 10000, alpha = 0.05, resampling = "paramBS", CPU,
                   seed, nested.levels.unique = FALSE, dec = 3){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  if(sum(grepl("cbind", formula)) != 0){
    stop("For data in wide format, please use function MANOVA.wide()")
  }
  
  test1 <- hasArg(CPU)
  if(!test1){
    CPU <- parallel::detectCores()
  }
  
  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }
  
  input_list <- list(formula = formula, data = data,
                     subject = subject, 
                     iter = iter, alpha = alpha, resampling = resampling, seed = seed)
  
  dat <- model.frame(formula, data)
  if (!(subject %in% names(data))){
    stop("The subject variable is not found!")
  }
  subject <- data[, subject]
  if (length(subject) != nrow(dat)){
    stop("There are missing values in the data.")
  }
  
  # no. of dimensions 
  p <- length(subject)/length(unique(subject))
  dat2 <- data.frame(dat, subject = subject)
  nf <- ncol(dat) - 1
  nadat <- names(dat)
  nadat2 <- nadat[-1]
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, (aa + 1)]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, (jj + 1)]))
  }
  lev_names <- expand.grid(levels)
  
  # number of hypotheses
  tmp <- 0
  for (i in 1:nf) {
    tmp <- c(tmp, choose(nf, i))
    nh <- sum(tmp)
  }
  
  if (nf == 1) {
    # one-way layout
    nest <- FALSE
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    dat2 <- dat2[order(dat2[, 2], dat2$subject), ]
    response <- dat2[, 1]    
    # contrast matrix
    hypo_matrices <- (diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(p)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(unique(subject)),
                     .drop = F)$Measure
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    MATS_out <- NA
    WTPS_out <- rep(NA, 2)
    quantiles <- matrix(NA, 2, 1)
    rownames(WTS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- MANOVA.Stat(data = response, n = n, hypo_matrices, iter = iter, alpha, resampling, n.groups = fl, p, CPU, seed, nf)    
    WTS_out[1, ] <- round(results$WTS, dec)
    MATS_out <- round(results$MATS, dec)
    WTPS_out <- round(results$WTPS, dec)
    quantiles <- results$quantiles
    names(quantiles) <- c("WTS_resampling", "MATS_resampling")
    mean_out <- matrix(round(results$Mean, dec), ncol = p, byrow = TRUE)
    Var_out <- results$Cov
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(nadat2, "n", rep("Means", p))   
    colnames(WTS_out) <- cbind ("Test statistic", "df",
                                "p-value")
    names(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
    #WTPS_out[WTPS_out == 0] <- "<0.001"
    colnames(MATS_out) <- "Test statistic"
    # end one-way layout ------------------------------------------------------
  } else {
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 2)]), ]
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    outcome_names <- rownames(nr_hypo)[1]  # names of outcome variables
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    gr <- nadat2[1]
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)
    n <- n$Measure/p
    
    # correct formula?
    if (length(fac_names) != nf && length(fac_names) != nh){
      stop("Something is wrong with the formula. Please specify all or no interactions.")
    }
    
    nested <- grepl(":", formula)
    nested2 <- grepl("%in%", formula)
    
    if (sum(nested) > 0 || sum(nested2) > 0){
      # nested
      nest <- TRUE
      
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
      nest <- FALSE
      
      ## adapting formula argument, if interaction term missing
      
      if (nrow(perm_names) != nh) {
        #stop("For crossed designs, an interaction term must be specified in the formula.")
        form2 <- as.formula(paste(outcome_names, "~", paste(fac_names, collapse = "*")))
        perm_names2 <- t(attr(terms(form2), "factors")[-1, ])
        fac_names2 <- attr(terms(form2), "term.labels")
        hyps <- HC_MANOVA(fl, perm_names2, fac_names2, p, nh)
        hypo_matrices <- hyps[[1]]
        fac_names2 <- hyps[[2]]
        # choose only relevant entries of the hypo matrices
        indices <- grep(":", fac_names2, invert = T)
        hypo_matrices <- lapply(indices, function(x) hypo_matrices[[x]])
        
      } else {
        hyps <- HC_MANOVA(fl, perm_names, fac_names, p, nh)
        hypo_matrices <- hyps[[1]]
        fac_names <- hyps[[2]]
      }
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
    
    
    n.groups <- prod(fl)
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
      results <- MANOVA.Stat(data = response, n, hypo_matrices[[i]],
                             iter, alpha, resampling, n.groups, p, CPU, seed, nf)
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
    colnames(descriptive) <- c(nadat2, "n", paste(rep("Mean", p), 1:p))
    colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
    colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
    # WTPS_out[WTPS_out == 0] <- "<0.001"
    colnames(MATS_out) <- "Test statistic"
    
  }
  # Output ------------------------------------------------------
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
  output$fl <- fl
  output$BSMeans <- results$BSmeans
  output$BSVar <- results$BSVar
  output$levels <- lev_names
  output$nested <- nest
  output$modelcall <- MANOVA
  
  # check for singular covariance matrix
  test <- try(solve(output$Covariance), silent = TRUE)
  if(!is.matrix(test)){
    warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
  }
  
  class(output) <- "MANOVA"
  return(output)
}
