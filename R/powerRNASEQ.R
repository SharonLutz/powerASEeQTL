#' @include helper_funcs.R
#' @importFrom MASS rnegbin glm.nb

powerRNASEQ <-
  function(n, mu, fold, phi, theta, n.simu, alpha, maf, methods, sim, propASE)
  {
    n2 = round(n*maf*maf) #number of subjects with both disease alleles
    n1 = round(n*2*maf*(1-maf)) #number of subjects with one disease allele
    n0 = n - n1 - n2 #number of subjects with no disease alleles
    x = rep(0:2, times=c(n0, n1, n2))
    
    pvalLin = pvalNB = pvalQP = pvalASE= pvalTreCASE = rep(1, n.simu)
    
    for(k in 1:n.simu){
      if(k%%10==0){print(paste("simulation",k,"of",n.simu))}
      
      # ------------------------------------------------
      # simulate data scenario 1: total read couont
      # ------------------------------------------------
      if(sim=="sim1"){
        y1 = y2 = rep(NA, n)
        
        y1 = rnegbin(n0, mu/2, phi) # X=0 A
        y1 = c(y1, rnegbin(n1, mu*(fold)/2, phi)) # X=1 B
        y1 = c(y1, rnegbin(n2, mu*(fold)/2, phi)) # X=2 B
        
        y2 = rnegbin(n0, mu/2, phi) # X=0 A
        y2 = c(y2, rnegbin(n1, mu/2, phi)) # X=1 A
        y2 = c(y2, rnegbin(n2, mu*(fold)/2, phi))# X=2 B
        
        y= y1+y2
        y1=round(y1*propASE)
        y2=round(y2*propASE)
        
        
      }#end of sim1
      
      # ------------------------------------------------
      # simulate data scenario 2: total read couont
      # ------------------------------------------------
      if(sim=="sim2"){
        
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
        
        nAS = round(propASE*y)
        
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
      }#end of sim2
      # ------------------------------------------------
      # different models for eQTL
      # ------------------------------------------------
      
      ## linear model
      yn = normscore(y)
      l1 = summary(lm(yn ~ x))
      pvalLin[k] = l1$coef[2,4]
      
      ## glm.nb
      g1 = glm.nb(y ~ x)
      s1 = summary(g1)
      pvalNB[k] = s1$coef[2,4]
      
      # Quasi Poisson 
      p1 = glm(y ~ x,family=quasipoisson(link=log))
      ps1 = summary(p1)
      pvalQP[k] = ps1$coef[2,4]
      
      
      # ------------------------------------------------
      # methods for ASE and both ASE eQTL
      # ------------------------------------------------
      if(sum(as.integer(methods=="eQTL.ASE"&methods=="ASE"))>0){
        nTotal = y1 + y2
        wkp    = which(nTotal >= 5)
        if(length(wkp) >= 5){
          ## trecase
          z2 = x
          z2[which(x==2)] = 4
          t2 = trecaseR(y, y1, y2, X=rep(1,n), z1=x, z2=z2)
          pvalTreCASE[k] = 1 - pchisq(t2$lrt,1)
          pvalASE[k] = 1 - pchisq(t2$lrtASE,1)
        }else{
          t1 = trecR(y, X=rep(1, n), z1=x, fam="negbin")
          pvalTreCASE[k]  = 1 - pchisq(t1$lrt,1)
          pvalASE[k]  = 1.0
        }
      }
      
      # ------------------------------------------------
      # methods both ASE eQTL only
      # ------------------------------------------------
      if(sum(as.integer(methods=="eQTL.ASE" & methods!="ASE"))>0){
        nTotal = y1 + y2
        wkp    = which(nTotal >= 5)
        if(length(wkp) >= 5){
          ## trecase
          z2 = x
          z2[which(x==2)] = 4
          t2 = trecaseR(y, y1, y2, X=rep(1,n), z1=x, z2=z2)
          pvalTreCASE[k] = 1 - pchisq(t2$lrt,1)
        }else{
          t1 = trecR(y, X=rep(1, n), z1=x, fam="negbin")
          pvalTreCASE[k]  = 1 - pchisq(t1$lrt,1)
        }
      }
      
      # ------------------------------------------------
      # methods ASE only
      # ------------------------------------------------
      if(sum(as.integer(methods!="eQTL.ASE" & methods=="ASE"))>0){
        nTotal = y1 + y2
        wkp    = which(nTotal >= 5)
        if(length(wkp) >= 5){
          ## trecase
          z2 = x
          z2[which(x==2)] = 4
          t2 = trecaseR(y, y1, y2, X=rep(1,n), z1=x, z2=z2)
          pvalASE[k] = 1 - pchisq(t2$lrtASE,1)
        }else{
          t1 = trecR(y, X=rep(1, n), z1=x, fam="negbin")
          pvalASE[k]  = 1.0
        }
      }
      
      
    }#end of simulation loop k
    # ------------------------------------------------
    
    ppLin = length(which(pvalLin < alpha))/n.simu
    ppNB = length(which(pvalNB < alpha))/n.simu
    ppQP = length(which(pvalQP < alpha))/n.simu
    
    if(sum(as.integer(methods=="ASE"))==1){
      ppASE = length(which(pvalASE < alpha))/n.simu
    }else{
      ppASE = NA
    }
    
    if(sum(as.integer(methods=="eQTL.ASE"))==1){
      ppTreCASE = length(which(pvalTreCASE < alpha))/n.simu
    }else{
      ppTreCASE = NA
    }

    # ------------------------------------------------
    
    # Return in this order: linear, negbin, quasi poisson, ase, trecase
    c(ppLin, ppNB, ppQP, ppASE, ppTreCASE)
    
  }

