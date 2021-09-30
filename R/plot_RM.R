# plotting routine for RM objects
# called by the utility functions

RMplotting <- function(plot.object, descr.object, gap, ax, col, pch,
                       legendpos, leg, leg.cex, ...){
  
  # extract relevant information first
  nf <- plot.object$nf
  lo <- grep("Lower", names(descr.object))
  up <- grep("Upper", names(descr.object))
  
  color <- col
  
  myplotci <- function(..., factor){
    plotrix::plotCI(...)
  }
  
  if(nf == 1){
    # one-way layouts
    levels <- plot.object$levels[[1]]
    y <- descr.object[, "Means"]
    li <- descr.object[, lo]
    ui <- descr.object[, up]
    xlab <- plot.object$fac_names
    
    myplotci(x = 1:length(levels), y, li = li,
             ui = ui, col = color[1], pch = pch[1], ...)
    if(ax){axis(side = 1, at = 1:1:length(levels), labels = levels)}
    
  } else if (nf == 2){
    # two-way plots
    levels <- plot.object$levels
    y <- matrix(descr.object[, "Means"], ncol = length(levels[[2]]), 
                  byrow = T)
    nlower <- matrix(descr.object[, lo], ncol = length(levels[[2]]), 
                     byrow = T)
    rownames(nlower) <- levels[[1]]
    colnames(nlower) <- levels[[2]]
    
    nupper <- matrix(descr.object[, up], ncol = length(levels[[2]]), 
                     byrow = T)
    rownames(nupper) <- levels[[1]]
    colnames(nupper) <- levels[[2]]

    myplotci(x = 1:length(levels[[2]]), y[1, ], li = nlower[1, ],
             ui = nupper[1, ],
             col = color[1], pch = pch[1],
             ...)
    if(ax){axis(side = 1, at = 1:1:length(levels[[2]]), labels = levels[[2]])}
    for(i in 2:nrow(y)){
      myplotci(x = (1:length(levels[[2]]))+gap*(i-1), y[i, ], li = nlower[i, ],
               ui = nupper[i, ], col = color[i], pch = pch[i], add = TRUE, ...)
    }
    if(leg){
      legend(legendpos, 
           legend = levels[[1]],
           col = color[1:nrow(y)],
           seg.len = 0.5, 
           lty = rep(1, nrow(y)), cex = leg.cex)
    }
  } else if (nf == 3){
    # three-way plots
    fl <- plot.object$fl
    levels <- plot.object$levels
    mu3 <- matrix(descr.object$Means, ncol = fl[3], byrow = TRUE)
    lower3 <- matrix(descr.object[, lo], ncol = fl[3], byrow = TRUE) 
    upper3 <- matrix(descr.object[, up], ncol = fl[3], byrow = TRUE)
    
    rownames(mu3) <- paste(rep(levels[[1]], each = length(levels[[2]])),
                           levels[[2]])
    colnames(mu3) <- levels[[3]]
    rownames(lower3) <- paste(rep(levels[[1]], each = length(levels[[2]])),
                              levels[[2]])
    colnames(lower3) <- levels[[3]]
    rownames(upper3) <- paste(rep(levels[[1]], each = length(levels[[2]])),
                              levels[[2]])
    colnames(upper3) <- levels[[3]]
  
    delta <- seq(from = 0, by = gap,
                 length = (fl[1] * fl[2] + 1))
    
    myplotci(x = 1:length(levels[[3]]), mu3[1, ], li = lower3[1, ],
             ui = upper3[1, ],
             col = color[1], pch = pch[1],
             ...)
    if(ax){axis(side = 1, at = 1:1:length(levels[[3]]), labels = levels[[3]])}
    
    color2 <- rep(color[1:fl[2]], fl[1])
    pch2 <- rep(pch[1:fl[1]], each = fl[2])
    
    for (j in 1:(fl[1]*fl[2]-1)){
      myplotci(x = ((1:fl[3]) + delta[j+1]),
               mu3[j+1, ], li = lower3[j+1, ],
               ui = upper3[j+1, ], add = TRUE, col = color2[j+1],
               pch = pch2[j+1], ...)
    }
    if(leg){
    legend(legendpos, 
           legend = c(plot.object$fac_names_original[1], levels[[1]], 
                      plot.object$fac_names_original[[2]], levels[[2]]),
           box.lty = 0,
           col = c(0, rep(1, fl[1]), 0, color[1:fl[2]]),
           pch = c(NA, pch[1:fl[1]], NA, rep(NA, fl[2])),
           seg.len = 0.5,
           lty = c(NA, rep(NA, fl[1]), NA, rep(1, fl[2])), cex = leg.cex)
    }
    
  } else if (nf >= 4) {
    stop("Higher-way interactions cannot be plotted!")
  }
}

