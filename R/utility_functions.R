xaxt <- NULL
#' Plot function for an RM object
#' 
#' Generic plot function for \code{RM} objects: Returns a plot of the mean 
#' values along with confidence intervals for a specified RM-model.
#' 
#' @param x An object of class \code{RM}
#' @param leg Logical: Should a legend be plotted?
#' @param ... Additional parameters to be passed to plot()
#' 

#' @export 
plot.RM <- function (x, leg = TRUE, ...) {
  
  object <- x
  dots <- list(...)
  a <- object$plotting
  b <- object$Descriptive
  nf <- a$nf
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
  
  label.x <- names(a$fl[nf])
  sum.fac <- sum(a$fl)
  xmax <- a$fl[nf]+ (sum.fac-a$fl[nf])*gap
  miny <- min(b$`Lower 95 % CI`)
  maxy <- max(b$`Upper 95 % CI`)
  # default values
  args <- list(a,
               b,
               lwd = 2, ylab = "Means", xlab = label.x,
               col = 1:(sum.fac-a$fl[nf]),
               pch = 1:18, legendpos = "topright",
               xlim = c(0.8, xmax + 0.1),
               ylim =c(miny, maxy), gap = 0.1, xaxt = "n", ax,
               leg, leg.cex = 0.5)
  
  args[names(dots)] <- dots
  do.call(RMplotting, args = args)
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
#' Returns a summary of the results including mean values, variances 
#' and sample sizes for all groups as well
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