#' Tests for Repeated Measures in Semi-Parametric Factorial Designs
#' 
#' The RM() function calculates the Wald-type statistic (WTS), the ANOVA-type 
#' statistic (ATS) as well as resampling versions of these test statistics for 
#' semi-parametric repeated measures designs.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variable and the right hand side contains the factor
#'   variables of interest. An interaction term must be specified. The time variable
#'   must be the last factor in the formula.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. Data must be in long format and must not contain missing values.
#' @param subject The column name of the subjects in the data. NOTE: Subjects within 
#' different groups of whole-plot factors must have individual labels, see Details for 
#' more explanation.
#' @param no.subf The number of sub-plot factors in the data, default is 1.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "Perm" (randomly permute 
#'    all observations), "paramBS" (parametric bootstrap approach) and "WildBS" 
#'    (wild bootstrap approach with Rademacher weights). Except for the Wild Bootstrap,
#'    all methods are applied to the WTS only.
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
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
#'  NOTE: The number of within-subject factors needs to be specified in the
#'  function call. If only one factor is present, it is assumed that this is a
#'  within-subjects factor (e.g. time).
#'  
#'  If subjects in different groups of the whole-plot factor have the same id, they will 
#'  not be identified as different subjects and thus it is erroneously assumed that their 
#'  measurements belong to one subject. Example: Consider a study with one whole-plot factor 
#'  "treatment" with leels verum and placebo aand one sub-plot factor "time" (4 measurements).
#'  If subjects in the placebo group are labelled 1-20 and subjects in the verum group have 
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
#'             subject = "Subject", no.subf = 2, iter = 1000, resampling = "Perm", CPU = 1)
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
#'                data = EEG, subject = "id", no.subf = 2, resampling = "WildBS",
#'                iter = 1000,  alpha = 0.01, CPU = 4, seed = 987, dec = 2)
#' summary(EEG_model)
#' }
#' 
#' @seealso \code{\link[GFD]{GFD}}, \code{\link[nparLD]{nparLD}}, \code{\link{MANOVA}}
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
               no.subf = 1, iter = 10000, alpha = 0.05, resampling = "Perm",
               CPU, seed, CI.method = "t-quantile", dec = 3){
  
   if (!(resampling %in% c("Perm", "paramBS", "WildBS"))){
     stop("Resampling must be one of 'Perm', 'paramBS' or 'WildBS'!")
   }
    
  if (!(CI.method %in% c("t-quantile", "resampling"))){
    stop("CI.method must be one of 't-quantile' or 'resampling'!")
  }
  
  input_list <- list(formula = formula, data = data,
                     subject = subject, 
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
  if (!(subject %in% names(data))){
    stop("The subject variable is not found!")
  }
  subject <- data[, subject]
  if (length(subject) != nrow(dat)){
    stop("There are missing values in the data.")
  }
  
  
  dat2 <- data.frame(dat, subject = subject)
  # check for missing values
  dt <- as.data.table(dat2)
  N <- NULL
  .N <- NULL
  if(NROW(dt[, .N, by = subject][unique(N)])!=1){
   stop("There are missing values in the data.")
  }
  
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
  
  # check that subjects are correctly labelled (at least for 1 sub-plot) factor
  if(no.subf == 1){
    if(nrow(data)/length(unique(subject)) != fl[length(fl)]){
      stop(paste0("The number of subjects (", length(unique(subject)), ") times the number of time points (", fl[length(fl)], ") does not equal the total number of observations (", nrow(data), ")."))
    }
  }

  if (nf == 1) {
    # one-way layout
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    # contrast matrix
    hypo <- diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    ATS_out <- matrix(NA, ncol = 4, nrow = 1)
    WTPS_out <- rep(NA, 2)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- RM.Stat.oneway(data = response, n = n, t = fl, hypo, iter = iter, 
                              alpha, resampling, seed, CI.method)
    WTS_out <- round(results$WTS, dec)
    ATS_out <- round(results$ATS, dec)
    WTPS_out <- round(results$WTPS, dec)
    mean_out <- round(results$Mean, dec)
    Var_out <- round(results$Cov, dec)
    CI <- round(results$CI, dec)
    colnames(CI) <- c("CIl", "CIu")
    descriptive <- cbind(lev_names, n, mean_out, CI)
    colnames(descriptive) <- c(nadat2, "n", "Means",
                               paste("Lower", 100 * (1 - alpha), "%", "CI"),
                               paste("Upper", 100 * (1 - alpha), "%", "CI"))
    
    names(WTS_out) <- cbind ("Test statistic", "df",
                                "p-value")
    names(ATS_out) <- cbind("Test statistic", "df1", "df2", "p-value")
    names(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"))
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$Covariance <- Var_out
    output$WTS <- WTS_out
    output$ATS <- ATS_out
    output$resampling <- WTPS_out
    output$plotting <- list(levels, fac_names, nf)
    names(output$plotting) <- c("levels", "fac_names", "nf")
    # end one-way layout ------------------------------------------------------
  } else {
    # no. of whole-plot (groups) and sub-plot (sub) factors
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 2)]), ]
    n.whole <- nf - no.subf
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    gr <- nadat2[1]
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)
    n.groups <- prod(fl[1:n.whole])
    n.sub <- prod(fl) / n.groups
    nind <- n[1, ]
    for (i in 1:(n.groups-1)){
      nind[i+1, ] <- n[(n.sub * i) + 1, ]
    }
    if (n.whole == 1){
    zzz <- unique(dat2[, gr])
    nind <- nind[order(order(zzz)), ]$Measure
    n <- rep(nind, each=n.sub)
    } else {
      nind <- nind$Measure
      n <- n$Measure
    }
    # no factor combinations with less than 2 observations
    if (0 %in% nind || 1 %in% nind) {
      stop("There is at least one factor-level combination
           with less than 2 observations!")
    }
    
    # number of whole-plot factors and their interactions
    tmp_whole <- 0
    for (i in 1:n.whole) {
      tmp_whole <- c(tmp_whole, choose(n.whole, i))
      whole_count <- sum(tmp_whole)
    }
    if (n.whole == 0){
      whole_count <- 0
    }
    
    # number of hypotheses
    tmp <- 0
    for (i in 1:nf) {
      tmp <- c(tmp, choose(nf, i))
      nh <- sum(tmp)
    }
    if (length(fac_names) != nh) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    
    hypo_matrices <- HC(fl, perm_names, fac_names)[[1]]
    fac_names <- HC(fl, perm_names, fac_names)[[2]]
    
    WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
    ATS_out <- matrix(NA, ncol = 4, nrow = length(hypo_matrices))
    WTPS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    rownames(WTPS_out) <- fac_names
    colnames(ATS_out) <- c("Test statistic", "df1", "df2", "p-value")
    # calculate results
    if (n.whole == 0 && nf !=1) {
      for (i in 1:length(hypo_matrices)) {
        results <- RM.Stat.sub(data = response, nind = nind[1], n, hypo_matrices[[i]],
                               iter, alpha, n.sub, n.groups, resampling, seed, CI.method)
        WTS_out[i, ] <- round(results$WTS, dec)
        ATS_out[i, ] <- round(results$ATS, dec)
        WTPS_out[i, ] <- round(results$WTPS, dec)
      }
    } else {  
      for (i in 1:length(hypo_matrices)) {
        results <- RM.Stat(data = response, nind, n, hypo_matrices[[i]],
                           iter, alpha, iii = i, whole_count, n.sub, n.groups, 
                           resampling, CPU, seed, CI.method)
        WTS_out[i, ] <- round(results$WTS, dec)
        ATS_out[i, ] <- round(results$ATS, dec)
        WTPS_out[i, ] <- round(results$WTPS, dec)
       # quant[i] <- results$quantile
      }
    }
    mean_out <- round(results$Mean, dec)
    Var_out <- round(results$Cov, dec)
    CI <- round(results$CI, dec)
    colnames(CI) <- c("CIl", "CIu")
    descriptive <- cbind(lev_names, n, mean_out, CI)
    colnames(descriptive) <- c(nadat2, "n", "Means",
                               paste("Lower", 100 * (1 - alpha),"%", "CI"),
                               paste("Upper", 100 * (1 - alpha),"%", "CI"))
    
    # calculate group means, variances and CIs ----------------------------
    mu <- list()
    sigma <- list()
    n_groups <- list()
    lower <- list()
    upper <- list()
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
    
    # Output ------------------------------------------------------
    colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
    colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"))
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$Covariance <- Var_out
    output$WTS <- WTS_out
    output$ATS <- ATS_out
    output$resampling <- WTPS_out
    output$plotting <- list(levels, fac_names, nf, no.subf, mu, lower, upper,
                            fac_names_original, dat2, fl, alpha, nadat2, lev_names)
    names(output$plotting) <- c("levels", "fac_names", "nf", "no.subf", "mu",
                                "lower", "upper", "fac_names_original", "dat2", "fl",
                                "alpha", "nadat2", "lev_names")
  }
  
  # check for singular covariance matrix
  test <- try(solve(output$Covariance), silent = TRUE)
  if(!is.matrix(test)){
    warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
  }
  
  class(output) <- "RM"
  return(output)
}
