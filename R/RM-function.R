#' Tests for Repeated Measures in Semi-Parametric Factorial Designs
#' 
#' The RM() function calculates the Wald-type statistic (WTS), the ANOVA-type 
#' statistic (ATS) as well as resampling versions of these test statistics for 
#' semi-parametric repeated measures designs.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variable and the right hand side contains the factor
#'   variables of interest. The within-subject factor(s)
#'   must be the last factor(s) in the formula, e.g. 
#'   \code{outcome ~ between1 * between2 * within1 * within2}.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. Data must be in long format and must not contain missing values.
#' @param subject The column name of the subjects in the data. NOTE: Subjects within 
#' different groups of between-subject factors must have individual labels, see Details for 
#' more explanation.
#' @param within Specifies the within-subject factor(s) in the formula. Either this
#'  or \code{no.subf} must be specified.
#' @param no.subf The number of within-subject factors in the data. Must be specified if
#' \code{within} is omitted.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "Perm" (randomly permute 
#'    all observations), "paramBS" (parametric bootstrap approach) and "WildBS" 
#'    (wild bootstrap approach with Rademacher weights). Except for the Wild Bootstrap,
#'    all methods are applied to the WTS only.
#' @param para If parallel computing should be used. Default is FALSE.
#' @param CPU The number of cores used for parallel computing. If not specified, cores
#'  are detected via \code{\link[parallel]{detectCores}}.
#' @param seed A random seed for the resampling procedure. If omitted, no 
#'   reproducible seed is set.
#' @param CI.method The method for calculating the quantiles used for the confidence intervals, 
#'  either "t-quantile" (the default) or "resampling" (the quantile of the resampled WTS).
#' @param dec Number of decimals the results should be rounded to. Default is 3.
#'   
#' @details The RM() function provides the Wald-type
#'  statistic as well as the ANOVA-type statistic for repeated measures designs
#'  with metric data as described in Friedrich et al. (2017).
#'  These are even applicable for non-normal error terms and/or heteroscedastic
#'  variances. It is implemented for designs with an arbitrary number of 
#'  between-subject (whole-plot) and within-subject (sub-plot) factors and 
#'  allows for different sample sizes. In addition to the
#'  asymptotic p-values, it also provides p-values based on resampling
#'  approaches.
#'  NOTE: The number of within-subject factors or their labels need
#'   to be specified in the function call. If only one factor is 
#'   present, it is assumed that this is a within-subject factor
#'   (e.g. time).
#'  
#'  If subjects in different groups of the between-subject factor have the same id, they will 
#'  not be identified as different subjects and thus it is erroneously assumed that their 
#'  measurements belong to one subject. Example: Consider a study with one between-subject factor 
#'  "treatment" with levels verum and placebo and one within-subject factor "time" (4 measurements).
#'  If subjects in the placebo group are labeled 1-20 and subjects in the verum group have 
#'  the same labels, the program erroneously assumes 20 individuals with 8 measurements each instead of 
#'  40 individuals with 4 measurements each.
#'   
#' @return An \code{RM} object containing the following components:
#' \item{Descriptive}{Some descriptive statistics of the data for all factor
#'   level combinations. Displayed are the number of individuals per factor
#'   level combination, the mean and 100*(1-alpha)\% confidence
#'   intervals (based on t-quantiles).}
#'  \item{Covariance}{The estimated covariance matrix.} 
#'  \item{WTS}{The value of the WTS along with degrees of freedom of the central 
#'  chi-square distribution and 
#'   corresponding p-value.}
#'  \item{ATS}{The value of the ATS, degrees of freedom of the central F distribution 
#'  and the corresponding p-value.}
#'  \item{resampling}{p-values for the test statistics based on the chosen resampling
#'   approach.}
#' 
#' @examples data(o2cons)
#' \dontrun{
#' oxy <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons, 
#'             subject = "Subject", no.subf = 2, iter = 1000,
#'             resampling = "Perm")
#' summary(oxy)
#' plot(oxy, factor = "Group") 
#'  
#' # For more details including the output of the examples also refer to the 
#' # package vignette.
#' 
#' # using the EEG data, consider additional within-subjects factors 'brain region' 
#' # and 'feature'
#' data(EEG)

#' EEG_model <- RM(resp ~ sex * diagnosis * feature * region, 
#'                data = EEG, subject = "id", within = c("feature", "region"),
#'                resampling = "WildBS",
#'                iter = 1000,  alpha = 0.01, seed = 987, dec = 2)
#' summary(EEG_model)
#' }
#' 
#' @seealso \code{\link[GFD]{GFD}}, \code{\link[nparLD]{nparLD}}, \code{\link{MANOVA}}
#' 
#' @references 
#' Friedrich, S., Konietschke, F., and Pauly, M. (2019). Resampling-Based Analysis
#'  of Multivariate Data and Repeated Measures Designs with the R Package MANOVA.RM.
#'  The R Journal, 11(2), 380-400.
#' 
#'  Friedrich, S., Brunner, E. and Pauly, M. (2017). Permuting longitudinal
#'  data in spite of the dependencies. Journal of Multivariate Analysis, 153, 255-265.
#' 
#'  Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and 
#'  Hoeller, Y. (2018). Testing Mean Differences among Groups: Multivariate and Repeated 
#'  Measures Analysis with Minimal Assumptions. Multivariate Behavioral Research, 53(3), 348-359,
#'  Doi: 10.1080/00273171.2018.1446320.
#'  
#'  Friedrich, S., Konietschke, F., Pauly, M. (2017). GFD - An 
#'  R-package for the Analysis of General Factorial Designs. 
#'  Journal of Statistical Software, 79(1), 1-18.
#' 
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

RM <- function(formula, data, subject,
               within, no.subf, iter = 10000, alpha = 0.05, resampling = "Perm",
               para = FALSE, CPU, seed, CI.method = "t-quantile", dec = 3){
  
  if (!(resampling %in% c("Perm", "paramBS", "WildBS"))){
    stop("Resampling must be one of 'Perm', 'paramBS' or 'WildBS'!")
  }
  
  if (!(CI.method %in% c("t-quantile", "resampling"))){
    stop("CI.method must be one of 't-quantile' or 'resampling'!")
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
                     subject = subject, 
                     iter = iter, alpha = alpha, resampling = resampling, seed = seed)
  #----------------------------------------------------------------------------------#
  # Determine names of within-subject factors if not given
  test3 <- hasArg(within)
  if(!test3){
    facs <- rownames(attr(terms(formula), "factors"))
    within <- facs[(length(facs)-no.subf+1):length(facs)]  
  }
  
  #----------------------------------------------------------------------------------#
  # prepare data for the calculations
  prepdt <- prepare.data(formula, data, subject, within)
  # extract relevant info
  EF <- prepdt[["EF"]]
  nf <- prepdt[["nf"]]
  fac_names_original <- prepdt[["fac_names_original"]]
  no.whole <- prepdt[["no.whole"]]
  lev.sub <- prepdt[["lev.sub"]]
  hypo_matrices <- prepdt[["hypo"]]
  fac_names <- prepdt[["fac"]]
  Yw2 <- prepdt[["data"]]
  nind <- prepdt[["n"]]
  fl <- prepdt[["fl"]]
  lev_names <- prepdt[["lev_names"]]
  n <- rep(nind, each = lev.sub)
  no.subf <- prepdt[["no.subf"]]
  # number of whole-plot factors and their interactions
  tmp_whole <- 0
  for (i in 1:no.whole) {
    tmp_whole <- c(tmp_whole, choose(no.whole, i))
    whole_count <- sum(tmp_whole)
  }
  if (no.whole == 0){
    whole_count <- 0
  }
  
  #---------------------------------------------------------------------#  
  WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
  ATS_out <- matrix(NA, ncol = 4, nrow = length(hypo_matrices))
  WTPS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2)
  rownames(WTS_out) <- fac_names
  rownames(ATS_out) <- fac_names
  rownames(WTPS_out) <- fac_names
  colnames(ATS_out) <- c("Test statistic", "df1", "df2", "p-value")
  colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
  colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"))
  
  # calculate results
  for (i in 1:length(hypo_matrices)) {
    results <- RM.Stat(Y = Yw2, nind, hypo_matrices[[i]],
                       iter, alpha, iii = i, whole_count, n.sub = lev.sub,
                       resampling, para, CPU, seed, CI.method)
    WTS_out[i, ] <- round(results$WTS, dec)
    ATS_out[i, ] <- round(results$ATS, dec)
    WTPS_out[i, ] <- round(results$WTPS, dec)
  }
  mean_out <- round(results$Mean, dec)
  Var_out <- round(results$Cov, dec)
  CI <- round(results$CI, dec)
  colnames(CI) <- c("CIl", "CIu")
  descriptive <- cbind(lev_names, n, mean_out, CI)
  colnames(descriptive) <- c(EF, "n", "Means",
                             paste("Lower", 100 * (1 - alpha),"%", "CI"),
                             paste("Upper", 100 * (1 - alpha),"%", "CI"))
  
  # calculate group means, variances and CIs ----------------------------
  mu <- list()
  sigma <- list()
  n_groups <- list()
  lower <- list()
  upper <- list()
  dat2 <- prepdt[["dat2"]]
  for (i in 1:nf) {
    mu[[i]] <- c(by(dat2[, 1], dat2[, i + 1], mean))
    sigma[[i]] <- c(by(dat2[, 1], dat2[, i + 1], var))
    n_groups[[i]] <- c(by(dat2[, 1], dat2[, i + 1], length))
    lower[[i]] <- mu[[i]] - sqrt(sigma[[i]] / n_groups[[i]]) *
      qt(1 - alpha / 2, df = n_groups[[i]])
    # HIER NUR t-QUANTIL
    upper[[i]] <- mu[[i]] + sqrt(sigma[[i]] / n_groups[[i]]) *
      qt(1 - alpha / 2, df = n_groups[[i]])
  }
  
  # for output: make sure the right factors are determined 
  # as within-subjects factors
  facs <- rownames(attr(terms(formula), "factors"))
  within_out <- facs[(length(facs)-no.subf+1):length(facs)] 
  
  # Output ------------------------------------------------------
  output <- list()
  output$input <- input_list
  output$Descriptive <- descriptive
  output$Covariance <- Var_out
  output$WTS <- WTS_out
  output$ATS <- ATS_out
  output$resampling <- WTPS_out
  output$plotting <- list(levels = prepdt[["levels"]], fac_names = fac_names,
                          nf = nf, no.subf = no.subf, mu = mu, 
                          lower = lower, upper = upper,
                          fac_names_original = fac_names_original, dat2 = dat2, 
                          fl = fl, alpha = alpha, nadat2 = EF, lev_names = lev_names)
  output$withinfactors <- within_out
  
  # check for singular covariance matrix
  test <- try(solve(output$Covariance), silent = TRUE)
  if(!is.matrix(test)){
    warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
  }
  
  class(output) <- "RM"
  return(output)
}
