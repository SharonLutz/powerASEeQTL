gradLogH1 <-
function(par, nA, nTotal, zeta){
  theta = par[1]
  pi1   = par[2]
  
  gradPi1 = gradTh = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    
    if(zeta[i]){
      gri     = gradBB(ni, ni0, pi1, theta)
      gradTh  = gradTh  + gri[1]
      gradPi1 = gradPi1 + gri[2]
    }else{
      gri     = gradBB(ni, ni0, 0.5, theta)
      gradTh  = gradTh + gri[1]
    }
    
  }
  
  c(gradTh, gradPi1)
}
