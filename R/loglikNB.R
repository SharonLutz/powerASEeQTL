loglikNB <-
function(phi, mu, y){
  logL = 0.0
  vphi = 1/phi
  lik0 = vphi*log(vphi) - lgamma(vphi)
  
  for(i in 1:length(y)){
    yi   = y[i]
    mui  = mu[i]
    if(yi==0){
      logL = logL + vphi*log(vphi) - vphi*log(vphi + mui)
    }else{
      logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui) 
      logL = logL - (vphi+yi)*log(vphi+mui) + lik0
    }
  }
  
  logL
}
