grad.theta <-
function(theta, nA, nTotal, zeta, pi)
{
  grad = 0.0
  
  for(i in 1:length(nA)){
    if(nA[i] > 0){
      ks = 0:(nA[i] - 1)
      
      if(zeta[i] > 0){
        grad = grad + sum(ks/(pi  + ks*theta)) 
      }else{
        grad = grad + sum(ks/(0.5 + ks*theta)) 
      }
    }
    
    if(nA[i] < nTotal[i]){
      ks = 0:(nTotal[i] - nA[i] - 1)
      
      if(zeta[i] > 0){
        grad = grad + sum(ks/(1 - pi + ks*theta))
      }else{
        grad = grad + sum(ks/(0.5 + ks*theta))
      }
    }
    
    ks = 1:(nTotal[i] - 1)
    grad = grad - sum(ks/(1 + ks*theta))
    
  }
  
  grad
}
