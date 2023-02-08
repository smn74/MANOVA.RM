#' Tests for Repeated Measures in Multivariate Semi-Parametric Factorial Designs
#' 
#' The multRM() function calculates the Wald-type statistic (WTS) and the modified ANOVA-type 
#' statistic (MATS) as well as resampling versions of these test statistics for multivariate
#' semi-parametric repeated measures designs.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the matrix of response variables (using cbind()) and the right hand side 
#'   contains the factor variables of interest. The within-subject factors must be specified
#'    last in the formula, e.g. \code{cbind(outcome1, outcome2) ~ between1 * between2 * within1 * within2}.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. Data must be in long format and must not contain missing values.
#' @param subject The column name of the subjects in the data. NOTE: Subjects within 
#' different groups of between-subject factors must have individual labels, see Details for 
#' more explanation.
#' @param within Specifies the within-subject factor(s) in the formula.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "paramBS" (parametric bootstrap
#'   approach) and "WildBS" (wild bootstrap approach with Rademacher weights).
#' @param para If parallel computing should be used. Default is FALSE.
#' @param CPU The number of cores used for parallel computing. If not specified, cores
#'  are detected via \code{\link[parallel]{detectCores}}.
#' @param seed A random seed for the resampling procedure. If omitted, no 
#'   reproducible seed is set.
#' @param dec Number of decimals the results should be rounded to. Default is 3.
#'   
#' @details The multRM() function provides the Wald-type
#'  statistic as well as the modified ANOVA-type statistic (Friedrich and Pauly, 2018) for repeated measures
#'  designs with multivariate metric outcomes.
#'  These methods are even applicable for non-normal error terms and/or heteroscedastic
#'  variances. Implemented are designs with an arbitrary number of 
#'  between-subject (whole-plot) and within-subject (sub-plot) factors and the methods
#'  allow for different sample sizes. In addition to the
#'  asymptotic p-values, p-values based on resampling
#'  approaches are provided.
#'  NOTE: The within-subject factors need to be specified in the
#'  function call (\code{within =}).
#'  
#'  If subjects in different groups of the between-subject factor have the same id, they will 
#'  not be identified as different subjects and thus it is erroneously assumed that their 
#'  measurements belong to one subject. See \code{\link{RM}} for more explanations and an example.
#'   
#' @return A \code{MANOVA} object containing the following components:
#'  \item{Descriptive}{Some descriptive statistics of the data for all factor 
#'   level combinations. Displayed are the number of individuals per factor 
#'   level combination and the vector of means (one column per dimension).}
#'   \item{Covariance}{The estimated covariance matrix.} 
#'   \item{WTS}{The value of the WTS along with degrees of freedom of the
#'   central chi-square distribution and p-value.} 
#'   \item{MATS}{The value of the MATS.} 
#'   \item{resampling}{p-values for the test statistic based on the
#'   chosen resampling approach.}
#'  
#' @examples
#' \dontrun{
#' data(EEG)
#' library(tidyr)
#' eeg <- spread(EEG, feature, resp)
#' fit <- multRM(cbind(brainrate, complexity) ~ sex * region, data = eeg, 
#'               subject = "id", within = "region")
#' summary(fit)
#' }
#' 
#' @seealso \code{\link{RM}}, \code{\link{MANOVA}}
#' 
#' @references Friedrich, S., Brunner, E. and Pauly, M. (2017). Permuting longitudinal
#'  data in spite of the dependencies. Journal of Multivariate Analysis, 153, 255-265.
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
#'    Friedrich, S., and Pauly, M. (2018). MATS: Inference for potentially singular and
#'   heteroscedastic MANOVA. Journal of Multivariate Analysis, 165, 166-179.
#'
#' @importFrom graphics axis legend par plot title
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var cov rbinom
#' @importFrom utils read.table
#' @importFrom methods hasArg
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster parSapply detectCores
#' @importFrom data.table as.data.table
#' 
#' @export

multRM <- function(formula, data, subject, within,
                   iter = 10000, alpha = 0.05, resampling = "paramBS",
                   para = FALSE, CPU, seed, dec = 3){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  if(para){
    test1 <- hasArg(CPU)
    if(!test1){
      CPU <- parallel::detectCores()
    }
  }
  
  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }
  
  input_list <- list(formula = formula, data = data,
                     iter = iter, alpha = alpha, resampling = resampling, seed = seed)
  #--------------------------------------------------------------------------------#
  # prepare data for the calculations
  prepdt <- prepare.data(formula, data, subject, within)
  # extract relevant info
  p <- prepdt[["p"]]
  EF <- prepdt[["EF"]]
  nf <- prepdt[["nf"]]
  split3 <- prepdt[["split3"]]
  no.whole <- prepdt[["no.whole"]]
  Yw2 <- prepdt[["data"]]
  hypo_matrices <- prepdt[["hypo"]]
  fac_names <- prepdt[["fac"]]
  nind <- prepdt[["n"]]
  fl <- prepdt[["fl"]]
  lev_names <- prepdt[["lev_names"]]
  no.subf <- prepdt[["no.subf"]]
  #-------------------------------------------------------------------------#
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
    results <- multRM.Statistic(Yw2, nind, hypo_matrices[[i]], iter, alpha, resampling, 
                                para, CPU,
                                seed, p, t=prod(fl[within]))
    WTS_out[i, ] <- round(results$WTS, dec)
    WTPS_out[i, ] <- round(results$WTPS, dec)
    MATS_out[i] <- round(results$MATS, dec)
    quantiles[i, ] <- results$quantiles
  }
  mean_out <- matrix(round(results$Mean, dec), ncol = p, byrow = TRUE)
  Var_out <- results$Cov
  descriptive <- cbind(unique(lev_names), rep(nind, each =prod(fl[within])) , mean_out)
  colnames(descriptive) <- c(EF, "n", split3)
  rownames(descriptive) <- NULL
  colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
  colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
  #WTPS_out[WTPS_out == 0] <- "<0.001"
  colnames(MATS_out) <- "Test statistic"
  
  # Output ------------------------------------------------------
  output <- list()
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
  #output$nested <- nest
  output$other <- list(no.subf = no.subf, no.whole = no.whole, p = p, within = within)
  output$nested <- FALSE
  output$modelcall <- multRM
  output$modeltype <- "multRM"
  
  # check for singular covariance matrix
  test <- try(solve(output$Covariance), silent = TRUE)
  if(!is.matrix(test)){
    warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
  }
  
  class(output) <- "MANOVA"
  return(output)
}
