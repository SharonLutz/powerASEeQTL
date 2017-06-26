logBB <-
function(ni, ni0, pi, theta){
  
  tmp0 = lchoose(ni, ni0)
  
  if(ni0 > 0){
    seq1 = seq(0, ni0-1, by=1)
    for(k in seq1){
      tmp0 = tmp0 + log(pi + k*theta)
    }
  }
  
  if(ni0 < ni){
    seq2 = seq(0, ni-ni0-1, by=1)
    for(k in seq2){
      tmp0 = tmp0 + log(1 - pi + k*theta)
    }
  }
  
  seq3 = 1:(ni-1)
  for(k in seq3){
    tmp0 = tmp0 - log(1 + k*theta)
  }
  
  tmp0
}
