loglikTheta <-
function(theta, pi, nA, nTotal, zeta){
  
  sumL = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    
    if(zeta[i]){
      sumL = sumL + logBB(ni, ni0, pi, theta)
    }else{
      sumL = sumL + logBB(ni, ni0, 0.5, theta)
    }
  }
  
  sumL
}
