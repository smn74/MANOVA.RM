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
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
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
#'               subject = "id", within = "region", iter = 200)
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
                   CPU, seed, dec = 3){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  if(!is.data.frame(data)){
    data <- as.data.frame(data)
  }
  
  output <- list()
  
  test1 <- hasArg(CPU)
  if(!test1){
    CPU <- parallel::detectCores()
  }
  
  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }
  
  input_list <- list(formula = formula, data = data,
                     iter = iter, alpha = alpha, resampling = resampling, seed = seed)
  
  dat <- model.frame(formula, data)
  if (!(subject %in% names(data))){
    stop("The subject variable is not found!")
  }
  subject <- data[, subject]
  if (length(subject) != nrow(dat)){
    stop("There are missing values in the data.")
  }
  
  dat <- data.frame(dat, subject = subject)
  nr_hypo <- attr(terms(formula), "factors")
  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  fac_names <- colnames(nr_hypo)
  fac_names_simple <- colnames(perm_names)
  
  if(!all(within %in% fac_names)){
    stop(paste0("The within-subjects factor ", within[which(!(within %in% fac_names_simple))], " is not part of the formula."))
  }
  
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
  whole <- fac_names_simple[which(!fac_names_simple %in% within)]
  lev.sub <- prod(fl[within])
  
  # correct formula?
  if (length(fac_names) != nf && length(fac_names) != nh){
    stop("Something is wrong with the formula. Please specify all or no interactions in crossed designs.")
  }
  
  # check that subjects are correctly labeled
    if(nrow(data)/length(unique(subject)) != prod(fl[within])){
      stop(paste0("The number of subjects (", length(unique(subject)), ") times the number of within-subject factor levels
                  (", prod(fl[within]), ") does not equal the total number of observations (", nrow(data), "). 
                  There are missing values in the data."))
    }
  
  if (nf == 1) {
    # one-way layout
    dat2 <- dat[order(dat[, 2]), ]
    dat2 <- dat2[order(dat2[, "subject"]), ]
    #fac.groups <- dat2[, 2]
    hypo <- list(diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl))
    Y <- list(dat2)
    lev.sub <- fl
    names(fl) <- within
    # end one-way layout ------------------------------------------------------
  } else {
    dat2 <- dat[do.call(order, dat[, c(2:(no.whole+1),      #order whole-plot factors
                                       ncol(dat),           # order by subject
                                       ((no.whole+2):(no.whole+2+no.subf-1)))]), ] # order sub-plot factors
   # dat2 <- dat2[order(dat2[, "subject"]), ]
    #fac.groups <- do.call(list, dat2[, 2:(nf+1)])
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    if(length(whole) ==0){
      Y <- list(dat2)
    } else {
       Y<- split(dat2, dat2[, whole], lex.order = TRUE)
    }
  }
  nind <- sapply(Y, nrow)/lev.sub
  
  ## adapting formula argument, if interaction term missing
  if (nrow(perm_names) != nh) {
    form2 <- as.formula(paste(outcome_names, "~", paste(fac_names, collapse = "*")))
    perm_names2 <- t(attr(terms(form2), "factors")[-1, ])
    fac_names2 <- attr(terms(form2), "term.labels")
    hyps <- HC(fl, perm_names2, fac_names2)
    hypo_matrices <- hyps[[1]]
    fac_names2 <- hyps[[2]]
    # choose only relevant entries of the hypo matrices
    indices <- grep(":", fac_names2, invert = TRUE)
    hypo <- lapply(indices, function(x) hypo_matrices[[x]])
    
  } else if(nf !=1){
    hyp <- HC(fl, perm_names, fac_names)
    hypo <- hyp[[1]]
    fac_names <- hyp[[2]]
  }

  hypo_matrices <- lapply(hypo, function(x) x %x% diag(p))
  # correcting for "empty" combinations (if no interaction specified)
  n.groups <- prod(fl[whole])
  if(nf != 1 & length(Y) != n.groups){
    index <- NULL
    for(i in 1:length(Y)){
      if(nrow(Y[[i]]) == 0){
        index <- c(index, i)
      }
    }
    Y <- Y[-index]
  }
  
  Ywhole <- lapply(Y, function(x) x$response)
  if (p==1){
    Ywhole <- lapply(Ywhole, function(x) as.matrix(x))
  }
  
  Yw2 <- lapply(Ywhole, function(x) matrix(t(x), nrow = nrow(x)/prod(fl[within]), ncol = p*prod(fl[within]), byrow=TRUE))
  
  # ---------------------- error detection ------------------------------------
  # no factor combinations with less than 2 observations
  if (0 %in% nind || 1 %in% nind) {
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
    results <- multRM.Statistic(Yw2, nind, hypo_matrices[[i]], iter, alpha, resampling, CPU, seed, p, t=prod(fl[within]))
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
  
  # check for singular covariance matrix
  test <- try(solve(output$Covariance), silent = TRUE)
  if(!is.matrix(test)){
    warning("The covariance matrix is singular. The WTS provides no valid test statistic!")
  }
  
  class(output) <- "MANOVA"
  return(output)
}
