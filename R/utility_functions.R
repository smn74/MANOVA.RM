xaxt <- NULL
#' Plot function for an RM object
#' 
#' Generic plot function for \code{RM} objects: Returns a plot of the mean values along with confidence 
#' intervals for a factor (combination) specified by the user.
#' 
#' @param x An object of class \code{RM}
#' @param CI.info If CI.info = TRUE, the mean values and confidence limits of the considered
#' contrast are printed.
#' @param ... Additional parameters to be passed to plot()
#' 
#' @details An additional argument \code{factor} can be used to specify the factor(s) used for plotting in two- and higher-way
#' layouts. See the examples for details.
#' 

#' @export 
plot.RM <- function (x, CI.info = FALSE, ...) {
  
  object <- x
  dots <- list(...)
  a <- object$plotting
  b <- object$Descriptive
  fac.names <- a$fac_names
  exist <- hasArg(factor) 
  
  if(length(fac.names) != 1){
    if(!exist){
      print("Please choose the factor you wish to plot (for interaction type something like group1:group2) and confirm by pressing 'Enter'")
      Faktor <- scan("", what="character")
      while(length(Faktor)==0){
        print("Please enter the name of the factor you wish to plot!")
        Faktor <- scan("", what="character")
      }
    } else {
      Faktor <- dots$factor
    }
  } else {
    Faktor <- fac.names
  }
  
  match.arg(Faktor, fac.names)
  h <- helper(a, b, Faktor)
  
  # to automatically create axis if not specified by user
   exist2 <- hasArg(xaxt)
   ax <- TRUE
   if(exist2){
     ax <- FALSE
   }
  
  if(!(hasArg(gap))){
    gap <- 0.1
  } else {
    gap <- dots$gap
  }
  
  xmax <- ncol(h$y)+ nrow(h$y)*gap
  # default values
  args <- list(h,
               lwd = 2, ylab = "Means", xlab = h$xlab, col = 1:length(h$legend),
               pch = 1:18, legendpos = "topright", xlim = c(0.8, xmax + 0.1),
               ylim =c(min(h$li), max(h$ui)), gap = 0.1, xaxt = "n", ax)
  
  args[names(dots)] <- dots
  do.call(newplotting, args = args)
  
  if (CI.info == TRUE){
    CI.out <- list()
    CI.out$mean <- h$y
    CI.out$lower <- h$li
    CI.out$upper <- h$ui
    return(CI.out)
  }
  
}

#' Display MANOVA object
#' 
#' Returns a short summary of the results (test statistics with p-values)
#' 
#' @param x A MANOVA object
#' @param ... Additional parameters (currently not used)
#' 
#' @export
print.MANOVA <- function(x, ...) {
  object <- x
  a <- object$input
  # avoid printing zeros
  WTS <- object$WTS
  WTS[WTS[, "p-value"] == 0, "p-value"] <- "<0.001"
  res <- object$resampling
  res[res == 0] <- "<0.001"
  
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(WTS)
  # cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  #  print(x$ATS)
  cat("\n", "modified ANOVA-Type Statistic (MATS):", "\n", sep = "")
  print(object$MATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  print(res)
}


#' Summarizing a MANOVA object
#' 
#' Returns a summary of the results including mean vectors and sample sizes for all groups as well
#' as test statistics with degrees of freedom and p-values
#' 
#' @param object A MANOVA object
#' @param ... Additional parameters (currently not used)
#' 
#' @export
summary.MANOVA <- function (object, ...) {
  a <- object$input
  b <- object$other
  # avoid printing zeros
  WTS <- object$WTS
  WTS[WTS[, "p-value"] == 0, "p-value"] <- "<0.001"
  res <- object$resampling
  res[res == 0] <- "<0.001"
  
  cat("Call:", "\n")
  print(a$formula)
  
  if(!is.null(b$within)){
  cat("A multivariate repeated measures analysis with ", b$no.subf, "within-subject factor(s) (", b$within, ")and ", b$no.whole,
      "between-subject factor(s).", "\n")
}
  
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(WTS)
  cat("\n", "modified ANOVA-Type Statistic (MATS):", "\n", sep = "")
  print(object$MATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  print(res)
}

#' Display an RM object
#' 
#' Returns a short summary of the results (test statistics with p-values)
#' 
#' @param x An RM object
#' @param ... Additional parameters (currently not used)
#' 
#' @export
print.RM <- function(x, ...) {
  object <- x
  a <- object$input
  # avoid printing zeros
  WTS <- object$WTS
  ATS <- object$ATS
  WTS[WTS[, "p-value"] == 0, "p-value"] <- "<0.001"
  ATS[ATS[, "p-value"] == 0, "p-value"] <- "<0.001"
  res <- object$resampling
  res[res == 0] <- "<0.001"
  
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(ATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  if(x$input$resampling == "Perm"){
    print(res[, 1, drop = FALSE])
  } else {
    print(res)
  }
}

#' Summarizing an RM object
#' 
#' Returns a summary of the results including mean values, variances and sample sizes for all groups as well
#' as test statistics with degrees of freedom and p-values
#' 
#' @param object An RM object
#' @param ... Additional parameters (currently not used)
#' 
#' @export
summary.RM <- function (object, ...) {
  a <- object$input
  b <- object$plotting
  # avoid printing zeros
  WTS <- object$WTS
  ATS <- object$ATS
  WTS[WTS[, "p-value"] == 0, "p-value"] <- "<0.001"
  ATS[ATS[, "p-value"] == 0, "p-value"] <- "<0.001"
  res <- object$resampling
  res[res == 0] <- "<0.001"
  
  cat("Call:", "\n")
  print(a$formula)
  cat("A repeated measures analysis with", b$no.subf, "within-subject factor(s) (", paste(object$withinfactors, collapse = ","), 
      ") and", b$nf-b$no.subf, "between-subject factor(s).", "\n")
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(ATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  if(a$resampling == "Perm"){
    print(res[, 1, drop = FALSE])
  } else {
    print(res)
  }
}


#' @export 
plot.MANOVA <- function(x, ...){
  stop("There is no plotting routine for MANOVA objects.")
}