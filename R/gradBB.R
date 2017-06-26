gradBB <-
function(ni, ni0, pi, theta){
  
  grTh = grPi = 0
  
  if(ni0 > 0){
    seq1 = seq(0, ni0-1, by=1)
    for(k in seq1){
      xx = 1/(pi + k*theta)
      grPi = grPi + xx
      grTh = grTh + k*xx
    }
  }
  
  if(ni0 < ni){
    seq2 = seq(0, ni-ni0-1, by=1)
    for(k in seq2){
      xx = 1/(1 - pi + k*theta)
      grPi = grPi - xx
      grTh = grTh + k*xx
    }
  }
  
  seq3 = 1:(ni-1)
  for(k in seq3){
    grTh = grTh - k/(1 + k*theta)
  }
  
  c(grTh, grPi)
}
