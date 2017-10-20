xaxt <- NULL
#' @export 
plot.RM <- function (x, ...) {
  
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
               lwd = 2, ylab = "Means", xlab = h$xlab, col = 1:length(fac.names),
               pch = 1:18, legendpos = "topright", xlim = c(0.8, xmax + 0.1),
               ylim =c(min(h$li), max(h$ui)), gap = 0.1, xaxt = "n", ax)
  
  args[names(dots)] <- dots
  do.call(newplotting, args = args)
}

#' @export
print.MANOVA <- function(x, ...) {
  a <- x$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(x$WTS)
  # cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  #  print(x$ATS)
  cat("\n", "modified ANOVA-Type Statistic (MATS):", "\n", sep = "")
  print(x$MATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  print(x$resampling)
}


#' @export
summary.MANOVA <- function (object, ...) {
  a <- object$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(object$WTS)
  cat("\n", "modified ANOVA-Type Statistic (MATS):", "\n", sep = "")
  print(object$MATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  print(object$resampling)
}

#' @export
print.RM <- function(x, ...) {
  a <- x$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(x$WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(x$ATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  print(x$resampling)
}

#' @export
summary.RM <- function (object, ...) {
  a <- object$input
  cat("Call:", "\n")
  print(a$formula)
  cat("\n", "Descriptive:", "\n", sep = "")
  print(object$Descriptive)
  cat("\n", "Wald-Type Statistic (WTS):", "\n", sep = "")
  print(object$WTS)
  cat("\n", "ANOVA-Type Statistic (ATS):", "\n", sep = "")
  print(object$ATS)
  cat("\n", "p-values resampling:", "\n", sep = "")
  print(object$resampling)
}

#' @export 
plot.MANOVA <- function(x, ...){
  stop("There is no plotting routine for MANOVA objects.")
}