gradLogH0 <-
function(par, nA, nTotal){
  theta = par[1]
  pi0   = 0.5
  
  gradTh = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    gri  = gradBB(ni, ni0, pi0, theta)
    
    gradTh = gradTh + gri[1]
  }
  
  gradTh
}
