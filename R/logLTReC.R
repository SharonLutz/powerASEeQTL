logLTReC <-
function(bxj, y, x, mu, b0, phi, fam){
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)
  
  mu1[ww2] = mu[ww2] * exp(bxj - b0)
  mu1[ww1] = mu[ww1] * (1 + exp(bxj))/(1 + exp(b0))
  
  if(fam=="negbin"){
    logL = loglikNB(phi, mu1, y)
  }else{
    logL = 0.0
    for(i in 1:length(y)){
      logL = logL + dpois(y[i], mu1[i], log = TRUE)  
    }
  }
  logL
}
