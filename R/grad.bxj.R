grad.bxj <-
function(bxj, y, x, nA, nTotal, zeta, mu, b0, phi, theta)
{
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)
  
  mu1[ww2] = exp(log(mu[ww2]) + bxj - b0)
  mu1[ww1] = exp(log(mu[ww1]) + log((1 + exp(bxj))/2) - log((1 + exp(b0))/2))
  
  dg.dmu = y/mu1 - (1 + phi*y)/(1 + phi*mu1)
  
  
  dmu.db = rep(0, length(mu1))
  dmu.db[ww1] = mu1[ww1]*exp(bxj)/(1 + exp(bxj))
  dmu.db[ww2] = mu1[ww2]
  
  w2use  = which(zeta > 0)
  
  pi = exp(bxj)/(1 + exp(bxj))
  
  dh.dpi = 0.0
  for(w1 in w2use){
    if(nA[w1] > 0){
      ks = 0:(nA[w1] - 1)
      dh.dpi = dh.dpi + sum(1/(pi + ks*theta)) 
    }
    
    if(nA[w1] < nTotal[w1]){
      ks = 0:(nTotal[w1] - nA[w1] - 1)
      dh.dpi = dh.dpi - sum(1/(1 - pi + ks*theta))
    }
  }
  
  dpi.db = exp(bxj)/((1 + exp(bxj))^2)
  
  grad = sum(dg.dmu*dmu.db) + dh.dpi*dpi.db
  grad
}
