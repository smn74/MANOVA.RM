#### function for generating hypotheses matrices in the MANOVA setting
# nf: number of factors involved
# fl: vector of factorlevels, i.e. (a, b, c, ....) of size nf
# p: Number of dimensions

#--------------------- crossed designs ------------------------------
# HC_MANOVA delegates to HC (same algorithm; p triggers square One() and %x% diag(p)).
HC_MANOVA <- function(fl, perm_names, names, p){
  HC(fl, perm_names, names, p = p)
}

# ----------------------- nested designs -------------------------------------
HN_MANOVA <- function(fl, p){
  nf <- length(fl)
  # centering matrix
  P <- function(x){
    P <- diag(x) - matrix(1 / x, ncol = x, nrow = x)
    return(P)
  }
  # scaled one-matricess
  One <- function(x){
    I <- matrix(1/x, ncol = x, nrow = x)
    return(I)
  }
  # function for calculating the kronecker product of several matrices
  kp <- function(A) {
    kp <- A[[1]]
    for (i in 2: length(A)) {
      kp <- kp %x% A[[i]]
    }
    return(kp)
  }
  hypo <- vector("list", nf)
  if (nf == 2) {
    hypo[[1]] <- P(fl[1]) %x% One(fl[2]) %x% diag(p)
    hypo[[2]] <- diag(fl[1]) %x% P(fl[2]) %x% diag(p)
  } else if (nf == 3) {
    hypo[[1]] <- P(fl[1]) %x% One(fl[2]) %x% One(fl[3]) %x% diag(p)
    hypo[[2]] <- diag(fl[1]) %x% P(fl[2]) %x% diag(fl[3]) %x% diag(p)
    hypo[[3]] <- diag(fl[1]) %x% diag(fl[2]) %x% P(fl[3]) %x% diag(p)
  } else {
    stop("Nested designs with four or more factors are not implemented!")
  }
  return(hypo)
}
