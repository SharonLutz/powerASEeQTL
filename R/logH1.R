logH1 <-
function(par, nA, nTotal, zeta){
  theta = par[1]
  pi1   = par[2]
  
  sumL  = 0
  
  for(i in 1:length(nA)){
    ni  = nTotal[i]
    ni0 = nA[i]
    
    if(zeta[i]){
      sumL = sumL + logBB(ni, ni0, pi1, theta)
    }else{
      sumL = sumL + logBB(ni, ni0, 0.5, theta)
    }
  }
  
  sumL
}
