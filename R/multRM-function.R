#' Tests for Repeated Measures in Multivariate Semi-Parametric Factorial Designs
#' 
#' The multRM() function calculates the Wald-type statistic (WTS) and the modified ANOVA-type 
#' statistic (MATS) as well as resampling versions of these test statistics for multivariate
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
#'  "treatment" with levels verum and placebo and one sub-plot factor "time" (4 measurements).
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

multRM <- function(formula, data, subject, within,
               iter = 10000, alpha = 0.05, resampling = "paramBS",
               CPU, seed, dec = 3){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  input_list <- list(formula = formula, data = data,
                     iter = iter, alpha = alpha, resampling = resampling)
  output <- list()
  
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
  
  
  dat <- data.frame(dat, subject = subject)
  # check for missing values
  # dt <- as.data.table(dat)
  # N <- NULL
  # .N <- NULL
  # if(NROW(dt[, .N, by = subject][unique(N)])!=1){
  #   stop("There are missing values in the data.")
  # }
  
  
  nr_hypo <- attr(terms(formula), "factors")
  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  fac_names <- colnames(nr_hypo)
  fac_names_simple <- colnames(perm_names)
  
  outcome_names <- rownames(nr_hypo)[1]  # names of outcome variables
  # extract names of outcome variables
  if (grepl("cbind", outcome_names)){
    split1 <- strsplit(outcome_names, "(", fixed = TRUE)[[1]][-1]
    split2 <- strsplit(split1, ")", fixed = TRUE)[[1]]
    split3 <- strsplit(split2, ",")[[1]]
  } else {
    split3 <- outcome_names
  }
  
  EF <- rownames(nr_hypo)[-1]  # names of influencing factors
  nf <- length(EF)
  names(dat) <- c("response", EF, "subject")
  #no. dimensions
  p <- ncol(as.matrix(dat$response))
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
  
  names(fl) <- fac_names_simple
  
  # determine within- and between-subject factors
  no.subf <- length(within)
  no.whole <- nf-no.subf
  whole <- fac_names_simple[fac_names_simple!= within]
  lev.sub <- prod(fl[within])
  
  # correct formula?
  if (length(fac_names) != nf && length(fac_names) != nh){
    stop("Something is wrong with the formula. Please specify all or no interactions in crossed designs.")
  }
  
  # if (nf == 1) {
  #   # one-way layout
  #   dat2 <- dat[order(dat[, 2]), ]
  #   fac.groups <- dat2[, 2]
  #   hypo_matrices <- list((diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(p))
  #   # end one-way layout ------------------------------------------------------
  # } else {
    dat2 <- dat[do.call(order, dat[, 2:(nf + 1)]), ]
    dat2 <- dat2[order(dat2[, "subject"]), ]
    fac.groups <- do.call(list, dat2[, 2:(nf+1)])
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
  #}
  Y.split.wholeplot <- split(dat2, dat2[, whole])#, lex.order = TRUE)
  nind <- sapply(Y.split.wholeplot, nrow)/lev.sub
  
    ## adapting formula argument, if interaction term missing
    # if (nrow(perm_names) != nh) {
    #   #stop("For crossed designs, an interaction term must be specified in the formula.")
    #   form2 <- as.formula(paste(outcome_names, "~", paste(fac_names, collapse = "*")))
    #   perm_names2 <- t(attr(terms(form2), "factors")[-1, ])
    #   fac_names2 <- attr(terms(form2), "term.labels")
    #   hyps <- HC_MANOVA(fl, perm_names2, fac_names2, p, nh)
    #   hypo_matrices <- hyps[[1]]
    #   fac_names2 <- hyps[[2]]
    #   # choose only relevant entries of the hypo matrices
    #   indices <- grep(":", fac_names2, invert = T)
    #   hypo_matrices <- lapply(indices, function(x) hypo_matrices[[x]])
    #   
    # } else if(nf !=1){
      hypo <- HC(fl, perm_names, fac_names)[[1]]
      fac_names <- HC(fl, perm_names, fac_names)[[2]]
      hypo_matrices <- lapply(hypo, function(x) x %x% diag(p))
      # hyps <- HC_multRM(fl, perm_names, fac_names, p, nh)
      # hypo_matrices <- hyps[[1]]
      # fac_names <- hyps[[2]]
    #}
  #}
  
  # correcting for "empty" combinations (if no interaction specified)
  n.groups <- prod(fl)
  # if(nf != 1 & length(Y) != n.groups){
  #   index <- NULL
  #   for(i in 1:length(Y)){
  #     if(nrow(Y[[i]]) == 0){
  #       index <- c(index, i)
  #     }
  #   }
  #   Y <- Y[-index]
  # }
  # Ygroups <- lapply(Y.split.groups, function(x) x$response)
  # if (p==1){
  #   Ygroups <- lapply(Ygroups, function(x) as.matrix(x))
  # }

  Ywhole <- lapply(Y.split.wholeplot, function(x) x$response)
  if (p==1){
    Ywhole <- lapply(Ywhole, function(x) as.matrix(x))
  }
 
  Yw2 <- lapply(Ywhole, function(x) matrix(t(x), nrow = nrow(x)/fl[within], ncol = p*fl[within], byrow=TRUE))
  
  
  # ---------------------- error detection ------------------------------------
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
    results <- multRM.Statistic(Yw2, nind, hypo_matrices[[i]], iter, alpha, resampling, CPU, seed, p, t=fl[within])
    WTS_out[i, ] <- round(results$WTS, dec)
    WTPS_out[i, ] <- round(results$WTPS, dec)
    MATS_out[i] <- round(results$MATS, dec)
    quantiles[i, ] <- results$quantiles
  }
  # time needed for resampling calculations
  time <- results$time
  mean_out <- matrix(round(results$Mean, dec), ncol = p, byrow = TRUE)
  Var_out <- results$Cov
  descriptive <- cbind(unique(lev_names), nind, mean_out)
  colnames(descriptive) <- c(EF, "n", split3)
  rownames(descriptive) <- NULL
  colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
  colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(MATS)"))
  #WTPS_out[WTPS_out == 0] <- "<0.001"
  colnames(MATS_out) <- "Test statistic"
  
  # Output ------------------------------------------------------
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
  
  
  # check for singular covariance matrix
  test <- try(solve(output$Covariance), silent = TRUE)
  if(!is.matrix(test)){
    warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
  }
  
  class(output) <- "MANOVA"
  return(output)
}
