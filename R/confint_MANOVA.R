#' The conf.reg() function calculates confidence regions for contrasts in multivariate factorial designs.
#' In the two-dimensional case, confidence ellipsoids can be plotted via the generic \code{plot()} function.
#' 
#' @param object A \code{MANOVA} object.
#' @param nullhypo In scenarios with more than one factor, the null hypothesis, i.e., 
#' the contrast of interest must be specified.
#'  
#' @return A \code{confreg} object containing the following components: 
#'   \item{center}{The center of the confidence ellipsoid.}
#'   \item{scale}{The scaling factors for the axis of the confidence ellipsoid calculated as \eqn{\sqrt{\lambda*c/N}}, where \eqn{\lambda} are the
#'   eigenvalues, c denotes the bootstrap quantile and N is the total sample size. See Friedrich and Pauly (2017) for details.}
#'   \item{eigenvectors}{The corresponding eigenvectors, which determine the axes of the ellipsoid.}
#' 
#' @examples data(EEG)
#' EEG_mod <- MANOVA(resp ~ sex * diagnosis, 
#'                     data = EEG, subject = "id", resampling = "paramBS", 
#'                     alpha = 0.05, iter = 100, CPU = 1)
#' conf.reg(EEG_mod, nullhypo = "sex")
#' 
#' 
#' @references Friedrich, S., and Pauly, M. (2017). MATS: Inference for potentially singular and 
#'   heteroscedastic MANOVA. arXiv preprint arXiv:1704.03731.
#'
#' @export
conf.reg <- function(object, nullhypo){
  
  confreg <- list()

  meanvec <- as.vector(t(object$Means))
  covmat <- object$Covariance
  n <- object$Descriptive$n
  N <- sum(n)
  nf <- object$nf
  factors <- object$factors
  p <- object$p
  
  if (nf == 1){
  c_star <- object$quantile[2]
  hypo <- p*object$H
  confreg$nullhypo <- factors
} else {
  index <- which(factors == nullhypo)
  c_star <- object$quantile[index, 2]
  hypo <- p*object$H[[index]]
  confreg$nullhypo <- nullhypo
}

hypo2 <- hypo[1:p, ]

center <- hypo2%*%meanvec
D <- diag(diag(covmat))

lambda <- eigen(hypo2%*%D%*%t(hypo2))$values
eig <- eigen(hypo2%*%D%*%t(hypo2))$vectors

confreg$center <- center
confreg$scale <- sqrt(lambda*c_star/sum(n))
confreg$eigenvectors <- eig

confreg$dim <- p

class(confreg) <- "confreg"
return(confreg)
}

# ----------------------------------------------------------------------------------------------#

#' @export
print.confreg <- function(x, ...){
  cat("Center:", "\n")
  print(x$center)
  cat("\n", "Scale:", "\n", sep = "")
  print(x$scale)
  cat("\n", "Eigenvectors:", "\n", sep = "")
  print(x$eigenvectors)
}


#' @export
plot.confreg <- function(x, ...){
  
  object <- x
  dots <- list(...)
  
  if(object$dim != 2){
    stop("Confidence ellipsoids can only be plotted in designs with 2 dimensions!")
  }
  
  eig <- object$eigenvectors
    
  elly <- ellipse::ellipse(t(eig), centre = object$center, scale = object$scale)
  
  args <- list(x = elly, type = "l",  
               #ylab = as.expression(bquote(paste(mu, "12-", mu, "22"))),
               main = paste("Confidence ellipsoid for factor", object$nullhypo))
  args[names(dots)] <- dots
  
  # plot for 2 groups, 2 dimensions
  do.call(plot, args = args)
  points(object$center[1], object$center[2], pch = 23, cex = 1.2, col = 1, bg = 1)
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)
  
}