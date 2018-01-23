# helper function to extract relevant information for plotCI()
helper <- function(plot.object, descr.object, factor, ...){
  
  h <- list()
  
  nf <- plot.object$nf
  lo <- grep("Lower", names(descr.object))
  up <- grep("Upper", names(descr.object))
  
  # one-way 
  if(nf == 1){
    h$levels <- plot.object$levels[[1]]
    h$y <- descr.object[, "Means"]
    h$li <- descr.object[, lo]
    h$ui <- descr.object[, up]
    h$xlab <- plot.object$fac_names
    h$code <- "main"
  } else {
    Faktor <- strsplit(factor, ":")[[1]]
    levels <- plot.object$levels
    nadat2 <- plot.object$nadat2
    mu <- plot.object$mu
    lower <- plot.object$lower
    upper <- plot.object$upper
    fl <- plot.object$fl
    fac_names_original <- plot.object$fac_names_original
    dat2 <- plot.object$dat2
    alpha <- plot.object$alpha
    
    # main effects
    if(factor %in% nadat2){
      index <- which(factor == nadat2)
      h$levels <- levels[[index]]
      h$y <- matrix(mu[[index]], ncol = fl[index])
      h$li <- lower[[index]]
      h$ui <- upper[[index]]
      h$xlab <- nadat2[index]
      h$code <- "main"
    }
    
    # two-way plots
    fac_names_twofold <- plot.object$fac_names_original[ - (1:nf)]
    fac_names_twofold <- fac_names_twofold[1:choose(nf, 2)]
    
    if (factor %in% fac_names_twofold) {
      fak1 <- Faktor[1]
      fak2 <- Faktor[2]
      posi <- which(fac_names_original[1:nf] == fak1)
      posi2 <- which(fac_names_original[1:nf] == fak2)
      
      nmu <- matrix(by(dat2[, 1], dat2[, c(fak1, fak2)], mean),
                    nrow = fl[posi])
      rownames(nmu) <- levels[[posi]]
      colnames(nmu) <- levels[[posi2]]
      nsigma <- matrix(by(dat2[, 1], dat2[, c(fak1, fak2)], var),
                       nrow = fl[posi])
      nn_groups <- matrix(by(dat2[, 1], dat2[, c(fak1, fak2)],
                             length), nrow = fl[posi])
      nlower <- nmu - sqrt(nsigma/ nn_groups) *
        qt(1 - alpha / 2, df = nn_groups)
      rownames(nlower) <- levels[[posi]]
      colnames(nlower) <- levels[[posi2]]
      
      nupper <- nmu + sqrt(nsigma / nn_groups) *
        qt(1 - alpha / 2, df = nn_groups)
      rownames(nupper) <- levels[[posi]]
      colnames(nupper) <- levels[[posi2]]
      
      # output of helper function
      h$levels <- levels[[posi2]]
      h$y <- nmu
      h$li <- nlower
      h$ui <- nupper
      h$xlab <- fak2
      h$legend <- levels[[posi]]
      h$code <- "2way"
    } else if (length(Faktor) == 3) {
      # three-way plots
      fak1 <- Faktor[1]
      fak2 <- Faktor[2]
      fak3 <- Faktor[3]
      
      posi1 <- which(fac_names_original[1:nf] == fak1)
      posi2 <- which(fac_names_original[1:nf] == fak2)
      posi3 <- which(fac_names_original[1:nf] == fak3)
      
      if (nf == 3){          
        mu3 <- matrix(descr.object$Means, ncol = fl[posi3], byrow = TRUE)
        lower3 <- matrix(descr.object[, lo], ncol = fl[posi3], byrow = TRUE) 
        upper3 <- matrix(descr.object[, up], ncol = fl[posi3], byrow = TRUE) 
      } else {
        # c(t(...)) to get correct order of factors, analogous to case above
        mu3 <- matrix(c(t(matrix(by(dat2[, 1], dat2[, c(fak1, fak3, fak2)], mean),
                          nrow = fl[posi1]))), ncol=fl[posi3], byrow = TRUE)
        nsigma <- matrix(c(t(matrix(by(dat2[, 1], dat2[, c(fak1, fak3, fak2)], var),
                             nrow = fl[posi1]))), ncol=fl[posi3], byrow = TRUE)
        nn_groups <- matrix(c(t(matrix(by(dat2[, 1], dat2[, c(fak1, fak3, fak2)],
                                   length), nrow = fl[posi1]))), ncol=fl[posi3], byrow = TRUE)
        lower3 <- mu3 - sqrt(nsigma/ nn_groups) *
          qt(1 - alpha / 2, df = nn_groups)
        upper3 <- mu3 + sqrt(nsigma / nn_groups) *
          qt(1 - alpha / 2, df = nn_groups)          
      }
      rownames(mu3) <- paste(rep(levels[[posi1]], each = length(levels[[posi2]])), levels[[posi2]])
      colnames(mu3) <- levels[[posi3]]
      rownames(lower3) <- paste(rep(levels[[posi1]], each = length(levels[[posi2]])), levels[[posi2]])
      colnames(lower3) <- levels[[posi3]]
      rownames(upper3) <- paste(rep(levels[[posi1]], each = length(levels[[posi2]])), levels[[posi2]])
      colnames(upper3) <- levels[[posi3]]
      
      # ouput
      h$levels <- levels
      h$y <- mu3   
      h$li <- lower3
      h$ui <- upper3
      h$xlab <- fak3
      h$code <- "3way"
      h$posi <- c(posi1, posi2, posi3)
      h$fl <- fl
      h$fac_names_original <- fac_names_original
      h$legend <- c(fac_names_original[posi1], levels[[posi1]], 
                    fac_names_original[[posi2]], levels[[posi2]])
      
    } else if (length(Faktor) >= 4) {
      stop("Higher-way interactions cannot be plotted!")
    }
  }
  return(h)
}