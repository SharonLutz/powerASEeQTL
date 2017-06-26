
powerRNAseq <-
function(n, mu, fold, phi, theta, n.simu, alpha, maf)
{
  
  library(MASS)
  
  n2 = round(n*maf*maf) #number of subjects with both disease alleles
  n1 = round(n*2*maf*(1-maf)) #number of subjects with one disease allele
  n0 = n - n1 - n2 #number of subjects with no disease alleles
  x = rep(0:2, times=c(n0, n1, n2))
  
  pval1 = pval2 = pval3 = pval4 = pval5 = pval6 = pval7 = useASE = rep(1, n.simu)

  for(k in 1:n.simu){
  if(k%%10==0){print(k)}
    
    # ------------------------------------------------
    # simulate data: total read couont
    # ------------------------------------------------
    
    b = log(fold)
    r = log(1+exp(b)) - log(2)
    y = rnegbin(n0, mu, phi)
    y = c(y, rnegbin(n1, mu*(1+fold)/2, phi))
    y = c(y, rnegbin(n2, mu*fold, phi))
    
    # ------------------------------------------------
    # simulate data: allele specific read
    # ------------------------------------------------

    pi0 = 0.5
    pi1 = exp(b)/(1 + exp(b))
    
    alpha0 = pi0/theta
    beta0  = 1/theta - alpha0
    
    alpha1 = pi1/theta
    beta1  = 1/theta - alpha1
    
    nAS = round(0.005*y)
    
    y1 = y2 = rep(NA, n)
    
    for(i in 1:n0){
      if(nAS[i] == 0){
        y1[i] = 0
        y2[i] = 0
      }else{
        y2[i] = rbetabinom.ab(1, nAS[i], alpha0, beta0)
        y1[i] = nAS[i] - y2[i]
      }
    }
    
    if(n1 > 0){
      for(i in (n0+1):(n0+n1)){
        if(nAS[i] == 0){
          y1[i] = 0
          y2[i] = 0
        }else{
          y2[i] = rbetabinom.ab(1, nAS[i], alpha1, beta1)
          y1[i] = nAS[i] - y2[i]
        }
      }
    }
    
    if(n2 > 0){
      for(i in (n0+n1+1):(n0+n1+n2)){
        if(nAS[i] == 0){
          y1[i] = 0
          y2[i] = 0
        }else{
          y2[i] = rbetabinom.ab(1, nAS[i], alpha0, beta0)
          y1[i] = nAS[i] - y2[i]
        }
      }    
    }
    
    # print(y1)
    # print(sum(y1))
    # print(y2)
    # print(sum(y2))
    
    # ------------------------------------------------
    # fit different models
    # ------------------------------------------------
    
    ## linear model
    yn = normscore(y)
    l1 = summary(lm(yn ~ x))
    pval1[k] = l1$coef[2,4]
    
    ## glm.nb
    g1 = glm.nb(y ~ x)
    s1 = summary(g1)
    pval2[k] = s1$coef[2,4]
    
    ## trec
    t1 = trecR(y, X=rep(1, n), z1=x, fam="negbin")
    pval3[k] = 1 - pchisq(t1$lrt,1)
    
    nTotal = y1 + y2
    wkp    = which(nTotal >= 5)
    if(length(wkp) >= 5){
      ## trecase
      z2 = x
      z2[which(x==2)] = 4
      t2 = trecaseR(y, y1, y2, X=rep(1,n), z1=x, z2=z2)
      pval4[k] = 1 - pchisq(t2$lrt,1)
      pval5[k] = 1 - pchisq(t2$lrtASE,1)
    }else{
      pval4[k]  = pval3[k]
      pval5[k]  = 1.0
      useASE[k] = 0
    }
  
    # ASEbinom pval6
    b1 <- binom.test(x = sum(y1), n = (sum(y1)+sum(y2)), p = 0.5, alternative = "two.sided", conf.level = (1-alpha))
    pval6[k] = b1$p.value
    
    # Poisson pval7
    p1 = glm(y ~ x,family=quasipoisson(link=log))
    ps1 = summary(p1)
    pval7[k] = ps1$coef[2,4]
    
  }
  
  pp1 = length(which(pval1 < alpha))/n.simu
  pp2 = length(which(pval2 < alpha))/n.simu
  pp3 = length(which(pval3 < alpha))/n.simu
  pp4 = length(which(pval4 < alpha))/n.simu
  pp5 = length(which(pval5 < alpha))/n.simu
  pp6 = length(which(pval6 < alpha))/n.simu
  pp7 = length(which(pval7 < alpha))/n.simu
  
  # pp1: linear model
  # pp2: glm.nb
  # pp3: TReC
  # pp4: TReCASE **I think
  # pp5: ASEchi2 **I think
  # pp6: ASEbinom
  # pp7: Poisson

  c(pp1, pp2, pp3, pp4, pp5, pp6, pp7)
}
