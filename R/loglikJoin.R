loglikJoin <-
function(bxj, y, x, nA, nTotal, zeta, mu, b0, phi, theta){
  
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)
  
  mu1[ww2] = exp(log(mu[ww2]) + bxj - b0)
  mu1[ww1] = exp(log(mu[ww1]) + log(1 + exp(bxj)) - log(1 + exp(b0)))
  
  pi1 = exp(bxj)/(1 + exp(bxj))
  par = c(theta, pi1)
  
  logTReC = loglikNB(phi, mu1, y)
  
  logASE  = logH1(par, nA, nTotal, zeta)
  
  logTReC + logASE
}
