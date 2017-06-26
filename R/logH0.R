logH0 <-
function(par, nA, nTotal){
  theta = par[1]
  pi0   = 0.5
  
  sumL = 0
  
  for(i in 1:length(nA)){
    ni   = nTotal[i]
    ni0  = nA[i]
    sumL = sumL + logBB(ni, ni0, pi0, theta)
  }
  
  sumL
}
