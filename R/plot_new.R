newplotting <- function(h, gap, ax, col, pch, legendpos, ...){
  
  color <- col
  
  myplotci <- function(..., factor){
    plotrix::plotCI(...)
  }
  
  if(h$code == "main"){
    # one-way and main effects
    myplotci(x = 1:length(h$levels), h$y, li = h$li,
             ui = h$ui, col = color[1], pch = pch[1], ...)
    if(ax){axis(side = 1, at = 1:1:length(h$levels), labels = h$levels)}
  } else if (h$code == "2way"){
    # two-way plots
    myplotci(x = 1:length(h$levels), h$y[1, ], li = h$li[1, ],
             ui = h$ui[1, ],
             col = color[1], pch = pch[1],
             ...)
    if(ax){axis(side = 1, at = 1:1:length(h$levels), labels = h$levels)}
    for(i in 2:nrow(h$y)){
      myplotci(x = (1:length(h$levels))+gap*(i-1), h$y[i, ], li = h$li[i, ],
               ui = h$ui[i, ], col = color[i], pch = pch[i], add = TRUE, ...)
    }
    legend(legendpos, 
           legend = h$legend,
           col = color[1:nrow(h$y)],
           seg.len = 0.5, 
           lty = rep(1, nrow(h$y)))
  } else if (h$code == "3way"){
    # three-way plots
    fl <- h$fl
    delta <- seq(from = 0, by = gap,
                 length = (fl[h$posi[1]] * fl[h$posi[2]] + 1))
    
    myplotci(x = 1:length(h$levels[[h$posi[3]]]), h$y[1, ], li = h$li[1, ],
             ui = h$ui[1, ],
             col = color[1], pch = pch[1],
             ...)
    if(ax){axis(side = 1, at = 1:1:length(h$levels[[h$posi[3]]]), labels = h$levels[[h$posi[3]]])}
    
    color2 <- rep(color[1:fl[h$posi[2]]], fl[h$posi[1]])
    pch2 <- rep(pch[1:fl[h$posi[1]]], each = fl[h$posi[2]])
    
    for (j in 1:(fl[h$posi[1]]*fl[h$posi[2]]-1)){
      myplotci(x = ((1:fl[h$posi[3]]) + delta[j+1]),
               h$y[j+1, ], li = h$li[j+1, ],
               ui = h$ui[j+1, ], add = TRUE, col = color2[j+1], pch = pch2[j+1], ...)
    }
    legend(legendpos, 
           legend = c(h$fac_names_original[h$posi[1]], h$levels[[h$posi[1]]], 
                      h$fac_names_original[[h$posi[2]]], h$levels[[h$posi[2]]]),
           box.lty = 0,
           col = c(0, rep(1, fl[h$posi[1]]), 0, color[1:fl[h$posi[2]]]),
           pch = c(NA, pch[1:fl[h$posi[1]]], NA, rep(NA, fl[h$posi[2]])),
           seg.len = 0.5,
           lty = c(NA, rep(NA, fl[h$posi[1]]), NA, rep(1, fl[h$posi[2]])))
    
  } 
}
