grad.bxj.trec <-
function(bxj, y, x, mu, b0, phi, fam)
{
  mu1 = mu
  ww2 = which(x==2)
  ww1 = which(x==1)

  mu1[ww2] = mu[ww2] * exp(bxj - b0)
  mu1[ww1] = mu[ww1] * (1 + exp(bxj))/(1 + exp(b0))
  
  dg.dmu = y/mu1 - (1 + phi*y)/(1 + phi*mu1)
  
  dmu.db = rep(0, length(mu1))
  dmu.db[ww1] = mu1[ww1]*exp(bxj)/(1 + exp(bxj))
  dmu.db[ww2] = mu1[ww2]
  
  grad = sum(dg.dmu*dmu.db)
  grad
}
