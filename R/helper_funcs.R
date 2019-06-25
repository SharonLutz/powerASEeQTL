
#' @title normscore
#' @param vec vec
normscore <-
  function(vec) {
    len  = length(na.omit(vec))+1
    rank = rank(na.omit(vec))
    ties = (rank - floor(rank)) > 0
    new.vec = vec[!is.na(vec)] 
    new.vec[!ties]=qnorm(rank[!ties]/len)
    new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
    vec[!is.na(vec)] = new.vec
    return(vec)
  }

#' @title grad.bxj.trec
#' @param bxj bxj
#' @param y y
#' @param x x
#' @param mu mu
#' @param b0 b0
#' @param phi phi
#' @param fam fam
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
    return(grad)
  }

#' @title loglikNB
#' @param phi phi
#' @param mu mu
#' @param y y
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
    
    return(logL)
  }

#' @title logLTReC
#' @param bxj bxj
#' @param y y
#' @param x x
#' @param mu mu
#' @param b0 b0
#' @param phi phi
#' @param fam fam
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
    return(logL)
  }

#' @title trecR
#' @param y y
#' @param x x
#' @param z1 z1
#' @param fam fam
#' @param nIter nIter
#' @param plotIt plotIt
#' @param trace trace
#' @param yfit yfit
trecR <-
  function(y, X, z1, fam, nIter=100, plotIt=FALSE, trace=FALSE, yfit=FALSE){
    
    #----------------------------------------------------------
    # initial model fitting
    #----------------------------------------------------------
    if(fam=="negbin"){
      g1      = glm.nb(y ~ X)
      phi0    = 1/g1$theta
      logLik0 = g1$twologlik
    }else{
      g1      = glm(y ~ X, family=poisson())
      phi0    = 0.0
      logLik0 = 2*logLik(g1)
    }
    
    b0   = 0.0
    mu0  = g1$fitted
    
    logLik = logLikOld =logLik0
    
    if(trace)
      message(sprintf("Initial: b0=%e, phi0=%e\n", b0, phi0))
    
    for(g in 1:nIter){
      
      parDiff = 0.0
      
      #--------------------------------------------------------
      # estimate bxj
      #--------------------------------------------------------
      
      if(plotIt){
        bs = seq(-0.3, 0.3, by=0.01)
        ll = rep(NA, length(bs))
        gs = rep(NA, length(bs))
        
        for(i in 1:length(bs)){
          bxj   = bs[i]
          gs[i] = grad.bxj.trec(bxj, y, x=z1, mu0, b0, phi0, fam)
          ll[i] = logLTReC(bxj, y, x=z1, mu0, b0, phi0, fam)
        }
        
        quartz()
        par(mfrow=c(2,1))
        plot(bs, gs, type="l")
        abline(h=0)
        plot(bs, ll, type="l")
      }
      
      o1 = optim(b0, fn=logLTReC, gr=grad.bxj.trec, y=y, x=z1, mu=mu0, b0=b0, phi=phi0, 
                 fam=fam, method="BFGS", control=list(fnscale=-1.0, trace=0))
      
      
      if(o1$convergence != 0){ 
        warning("g =", g, " fail BFGS\n")
        return(NULL)
      }
      
      if(plotIt){
        abline(v=o1$par)
      }
      
      b1 = o1$par
      grad.bxj.trec(b1, y, z1, mu0, b0, phi0, fam)
      
      mu1  = mu0
      ww2  = which(z1==2)
      ww1  = which(z1==1)
      
      mu1[ww2] = mu0[ww2] * exp(b1 - b0)
      mu1[ww1] = mu0[ww1] * (1 + exp(b1))/(1 + exp(b0))
      
      logLik1 = 2.0*o1$value
      
      if(trace){
        message(sprintf("\ng=%d", g))
        message(sprintf("  Estimate bxj: logLik_TReC=(%e, %e)", logLikOld, logLik1))
      }
      
      if((logLik1 - logLikOld)/abs(logLikOld) < -1e-4){stop("liklihood decreases for bxj\n")}
      
      if(parDiff < abs(b1 - b0)) parDiff = abs(b1 - b0)
      
      b0  = b1
      mu0 = mu1
      logLikOld = logLik1
      
      #-----------------------------------------------------------
      # estimate mu and phi
      #-----------------------------------------------------------
      
      bx1  = rep(0, length(y))
      ww1  = which(z1==1)
      ww2  = which(z1==2)
      bx1[ww1] = log((1 + exp(b0))/2)
      bx1[ww2] = b0
      
      if(fam=="negbin"){
        g1   = glm.nb(y ~ X + offset(bx1))
        phi1 = 1/g1$theta
        logLik1 = g1$twologlik
      }else{
        g1 = glm(y ~ X, family=poisson(), offset=bx1)
        phi1 = 0.0
        logLik1 = 2*logLik(g1)
      }
      
      mu1  = g1$fitted
      
      if((logLik1 - logLikOld)/abs(logLikOld) < -1e-5){stop("liklihood decreases for phi\n")}
      
      if(parDiff < abs(phi1 - phi0)) parDiff = abs(phi1 - phi0)
      
      if(trace){
        message(sprintf("  Estimate mu & phi: phi=%.2f, logLik_TReC=(%e, %e)", phi1, logLikOld, logLik1))
      }
      
      phi0 = phi1
      mu0  = mu1
      logLikOld = logLik1
      
      #-----------------------------------------------------------
      # check convergence
      #-----------------------------------------------------------
      
      logLik = c(logLik, logLik1)
      
      if(trace){
        message(sprintf("g=%d, parDiff=%.e, b0=%e, phi0=%e", 
                        g, parDiff, b0, phi0))
      }
      
      if(parDiff < 5e-5) break
      
    }
    
    if(yfit){
      l1 = list(b=b0, phi=phi0, logLik=logLik, lrt=logLik1 - logLik0, fitted=g1$fitted)
    }else{
      l1 = list(b=b0, phi=phi0, logLik=logLik, lrt=logLik1 - logLik0)
    }
    
    return(l1)
  }

#' @title logBB
#' @param ni ni
#' @param ni0 ni0
#' @param pi pi
#' @param theta theta
logBB <-
  function(ni, ni0, pi, theta){
    
    tmp0 = lchoose(ni, ni0)
    
    if(ni0 > 0){
      seq1 = seq(0, ni0-1, by=1)
      for(k in seq1){
        tmp0 = tmp0 + log(pi + k*theta)
      }
    }
    
    if(ni0 < ni){
      seq2 = seq(0, ni-ni0-1, by=1)
      for(k in seq2){
        tmp0 = tmp0 + log(1 - pi + k*theta)
      }
    }
    
    seq3 = 1:(ni-1)
    for(k in seq3){
      tmp0 = tmp0 - log(1 + k*theta)
    }
    
    return(tmp0)
  }

#' @title logH0
#' @param par par
#' @param nA nA
#' @param nTotal nTotal
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
    
    return(sumL)
  }

#' @title gradBB
#' @param ni ni
#' @param ni0 ni0
#' @param pi pi
#' @param theta theta
gradBB <-
  function(ni, ni0, pi, theta){
    
    grTh = grPi = 0
    
    if(ni0 > 0){
      seq1 = seq(0, ni0-1, by=1)
      for(k in seq1){
        xx = 1/(pi + k*theta)
        grPi = grPi + xx
        grTh = grTh + k*xx
      }
    }
    
    if(ni0 < ni){
      seq2 = seq(0, ni-ni0-1, by=1)
      for(k in seq2){
        xx = 1/(1 - pi + k*theta)
        grPi = grPi - xx
        grTh = grTh + k*xx
      }
    }
    
    seq3 = 1:(ni-1)
    for(k in seq3){
      grTh = grTh - k/(1 + k*theta)
    }
    
    return(c(grTh, grPi))
  }


#' @title gradLogH0
#' @param par par
#' @param nA nA
#' @param nTotal nTotal
gradLogH0 <-
  function(par, nA, nTotal){
    theta = par[1]
    pi0   = 0.5
    
    gradTh = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      gri  = gradBB(ni, ni0, pi0, theta)
      
      gradTh = gradTh + gri[1]
    }
    
    return(gradTh)
  }

#' @title logH1
#' @param par par
#' @param nA nA
#' @param nTotal nTotal
#' @param zeta zeta
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
    
    return(sumL)
  }

#' @title gradLogH1
#' @param par par
#' @param nA nA
#' @param nTotal nTotal
#' @param zeta zeta
gradLogH1 <-
  function(par, nA, nTotal, zeta){
    theta = par[1]
    pi1   = par[2]
    
    gradPi1 = gradTh = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      
      if(zeta[i]){
        gri     = gradBB(ni, ni0, pi1, theta)
        gradTh  = gradTh  + gri[1]
        gradPi1 = gradPi1 + gri[2]
      }else{
        gri     = gradBB(ni, ni0, 0.5, theta)
        gradTh  = gradTh + gri[1]
      }
      
    }
    
    return(c(gradTh, gradPi1))
  }


#' @title aseR
#' @param nA nA
#' @param nTotal nTotal
#' @param zeta zeta
#' @param maxIt maxIt
#' @param trace trace
aseR <-
  function(nA, nTotal, zeta, maxIt=50, trace=0) {
    
    # --------------------------------------------------------------    
    # check input
    # --------------------------------------------------------------
    
    if(!is.numeric(nA)){
      stop("nA must be a numerical vector\n")
    }
    
    if(!is.numeric(nTotal)){
      stop("nTotal must be a numerical vector\n")
    }
    
    if(length(nA) != length(nTotal)){
      stop("nA and nTotal have different lengths\n")  
    }
    
    if(any(nA > nTotal)){
      stop("nA must be smaller than nTotal\n")  
    }
    
    # -------------------------------------------------------------
    # N is the number of samples
    # -------------------------------------------------------------
    
    N = length(nA)
    
    # -------------------------------------------------------------
    # First, find MLE under H0
    # -------------------------------------------------------------
    
    par0 = 0.1
    op0  = optim(par0, logH0, gr=gradLogH0, nA=nA, nTotal=nTotal, 
                 method="L-BFGS-B", lower=1e-16, upper=Inf, 
                 control=list(fnscale=-1))
    
    theta0 = op0$par
    
    # -------------------------------------------------------------
    # Next, find MLE under H1
    # -------------------------------------------------------------
    
    par1 = c(theta0, 0.5)
    op1  = optim(par1, logH1, gr=gradLogH1, nA=nA, nTotal=nTotal, zeta=zeta, 
                 method="L-BFGS-B", lower=c(0,0) + 1e-16, upper=c(Inf, 1-1e-16), 
                 control=list(fnscale=-1))
    
    return(list(parH0=theta0,  logLikH0=op0$value, 
         parH1=op1$par, logLikH1=op1$value))
  }

#' @title grad.bxj
#' @param bxj bxj
#' @param y y
#' @param x x
#' @param nA, nA
#' @param nTotal nTotal
#' @param zeta zeta
#' @param mu mu
#' @param b0 b0
#' @param phi phi
#' @param theta theta
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
    return(grad)
  }

#' @title loglikJoin
#' @param bxj bxj
#' @param y y
#' @param x x
#' @param nA, nA
#' @param nTotal nTotal
#' @param zeta zeta
#' @param mu mu
#' @param b0 b0
#' @param phi phi
#' @param theta theta
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
    
    return(logTReC + logASE)
  }

#' @title loglikTheta
#' @param theta theta
#' @param pi pi
#' @param nA, nA
#' @param nTotal nTotal
#' @param zeta zeta
loglikTheta <-
  function(theta, pi, nA, nTotal, zeta){
    
    sumL = 0
    
    for(i in 1:length(nA)){
      ni   = nTotal[i]
      ni0  = nA[i]
      
      if(zeta[i]){
        sumL = sumL + logBB(ni, ni0, pi, theta)
      }else{
        sumL = sumL + logBB(ni, ni0, 0.5, theta)
      }
    }
    
    return(sumL)
  }

#' @title do_ASE
#' @param y y
#' @param y1 y1
#' @param y2 y2
#' @param x x
#' @param z1 z1
#' @param z2 z2
do_ASE <- function(y, y1, y2, X, z1, z2){
  #----------------------------------------------------------
  # initial model fitting
  #---------------------------------------------------------- 
  
  g0 = glm.nb(y ~ X)
  g0
  
  nTotal = as.numeric(y1 + y2)
  nA     = as.numeric(y2)
  nA[z2==3] = y1[z2==3]
  
  wkp    = which(nTotal >=5)
  if(length(wkp) < 5){
    stop("no enough allele specific reads\n")  
  }
  
  nTotal = nTotal[wkp]
  nA     = nA[wkp]
  
  wkp1   = which(abs(z2[wkp] - 2) < 1.5)      
  zeta   = rep(0, length(nTotal))
  zeta[wkp1] = 1
  
  a1 = aseR(nA, nTotal, zeta)
  return(2.0*(a1$logLikH1 - a1$logLikH0))
}

#' @title trecaseR
#' @param y y
#' @param y1 y1
#' @param y2 y2
#' @param x x
#' @param z1 z1
#' @param z2 z2
#' @param plotIt plotIt
#' @param traceIt traceIt
trecaseR <-
  function(y, y1, y2, X, z1, z2, plotIt=FALSE, traceIt=FALSE){
    
    #----------------------------------------------------------
    # initial model fitting
    #---------------------------------------------------------- 
    
    g0 = glm.nb(y ~ X)
    g0
    
    nTotal = as.numeric(y1 + y2)
    nA     = as.numeric(y2)
    nA[z2==3] = y1[z2==3]
    
    wkp    = which(nTotal >=5)
    if(length(wkp) < 5){
      stop("no enough allele specific reads\n")  
    }
    
    nTotal = nTotal[wkp]
    nA     = nA[wkp]
    
    wkp1   = which(abs(z2[wkp] - 2) < 1.5)      
    zeta   = rep(0, length(nTotal))
    zeta[wkp1] = 1
    
    a1 = aseR(nA, nTotal, zeta)
    a1
    
    pi0    = 0.5
    theta0 = a1$parH0
    phi0   = 1/g0$theta
    b0     = 0.0
    mu0    = g0$fitted
    b0_ase = log(pi0/(1-pi0))
    g1     = g0
    
    logLik = NULL
    
    if(traceIt){
      message(sprintf("Initial: b0=%e, pi0=%e, b0_ase=%e, theta0=%e, phi0=%e\n", 
                      b0, pi0, b0_ase, theta0, phi0))
    }
    
    par0       = c(theta0, pi0)
    logLikNULL = 2*logH1(par0, nA, nTotal, zeta) + 2*loglikNB(phi0, mu0, y)
    
    for(g in 1:100){
      parDiff = 0.0
      
      #-----------------------------------------------------------
      # estimate bxj
      #-----------------------------------------------------------
      
      if(plotIt){
        bs = seq(-1.0, 1.0, by=0.1)
        ll = rep(NA, length(bs))
        
        gs = rep(NA, length(bs))
        
        for(i in 1:length(bs)){
          bxj   = bs[i]
          gs[i] = grad.bxj(bxj, y, x=z1, nA, nTotal, zeta, mu0, b0, phi0, theta0)
          ll[i] = loglikJoin(bxj, y, x=z1, nA, nTotal, zeta, mu0, b0, phi0, theta0)
        }
        
        quartz()
        par(mfrow=c(2,1))
        plot(bs, gs, type="l")
        abline(h=0)
        plot(bs, ll, type="l")
      }
      
      par0    = c(theta0, pi0)
      logLik0 = 2*logH1(par0, nA, nTotal, zeta) + g1$twologlik
      
      op = optim(par=b0, fn=loglikJoin, gr=grad.bxj, y=y, x=z1, nA=nA, nTotal=nTotal, 
                 zeta=zeta, mu=mu0, b0=b0, phi=phi0, theta=theta0, method = "BFGS", 
                 control=list(fnscale=-1.0))
      
      b1  = op$par
      
      pi1 = exp(b1)/(1 + exp(b1))
      
      par1 = c(theta0, pi1)
      mu1  = mu0
      ww2  = which(z1==2)
      ww1  = which(z1==1)
      
      mu1[ww2] = mu0[ww2] * exp(b1 - b0)
      mu1[ww1] = mu0[ww1] * (1 + exp(b1))/(1 + exp(b0))
      
      logLik1 = 2*logH1(par1, nA, nTotal, zeta) + 2*loglikNB(phi0, mu1, y)
      
      if(traceIt){
        message(sprintf("\ng=%d", g))
        message(sprintf("  logLik_TReC=(%e, %e)", g1$twologlik, 2*loglikNB(phi0, mu1, y)))
        tmp0 = 2*logH1(par0, nA, nTotal, zeta)
        tmp1 = 2*logH1(par1, nA, nTotal, zeta)
        message(sprintf("  logLik_ASE =(%e, %e)", tmp0, tmp1))
      }
      
      if((logLik1 - logLik0)/abs(logLik0) < -0.01){stop("liklihood decreases for bxj\n")}
      
      if(parDiff < abs(b1 - b0)) parDiff = abs(b1 - b0)
      
      b0  = b1
      pi0 = pi1
      mu0 = mu1
      
      #-----------------------------------------------------------
      # estimate theta
      #-----------------------------------------------------------
      if(plotIt){
        ths = seq(0,0.1,by=0.001)
        gs  = rep(NA, length(ths))
        ll  = rep(NA, length(ths))
        
        for(i in 1:length(ths)){
          gs[i] = grad.theta(ths[i], nA, nTotal, zeta, pi0)
          ll[i] = logH1(par=c(ths[i], pi0), nA, nTotal, zeta)
        }
        
        quartz()
        par(mfrow=c(2,1))
        plot(ths, gs, type="l")
        abline(h=0)
        plot(ths, ll, type="l")
      }
      
      # uk = uniroot(grad.theta, c(0, 1), nA=nA, nTotal=nTotal, 
      #               zeta=zeta, pi=pi0, tol=1e-10)
      # theta1 = uk$root
      op = optimize(loglikTheta, interval=c(0, 1000), pi=pi0, nA=nA, 
                    nTotal=nTotal, zeta=zeta, maximum=TRUE)
      
      theta1 = op$maximum
      
      llase0 = 2*loglikTheta(theta0, pi0, nA, nTotal, zeta)
      llase1 = 2*loglikTheta(theta1, pi0, nA, nTotal, zeta)
      
      if(traceIt){
        message(sprintf("\ng=%d, estimate theta", g))
        message(sprintf("  (llase0, llase1)=(%e, %e)", llase0, llase1))
      }
      
      if((llase1 - llase0)/abs(llase0) < -0.01){stop("liklihood decreases for theta\n")}
      
      if(parDiff < abs(theta1 - theta0)) parDiff = abs(theta1 - theta0)
      theta0 = theta1
      
      #-----------------------------------------------------------
      # estimate mu and phi
      #-----------------------------------------------------------
      
      bx1  = rep(0, length(y))
      ww1  = which(z1==1)
      ww2  = which(z1==2)
      bx1[ww1] = log((1 + exp(b0))/2)
      bx1[ww2] = b0
      g1   = glm.nb(y ~ X + offset(bx1))
      mu1  = g1$fitted
      phi1 = 1/g1$theta
      
      lltrec0 = 2*loglikNB(phi0, mu0, y)
      lltrec1 = 2*loglikNB(phi1, mu1, y)
      
      if(traceIt){
        message(sprintf("\ng=%d, estimate phi", g))
        message(sprintf("  (lltrec0, lltrec1)=(%e, %e)", lltrec0, lltrec1))
      }
      
      if((lltrec1 - lltrec0)/abs(lltrec0) < -0.01){stop("liklihood decreases for phi")}
      
      if(parDiff < abs(phi1 - phi0)) parDiff = abs(phi1 - phi0)
      
      phi0 = phi1
      mu0  = mu1
      
      #-----------------------------------------------------------
      # check convergence
      #-----------------------------------------------------------
      
      par0 = c(theta0, pi0)
      logLik = c(logLik, 2*logH1(par0, nA, nTotal, zeta) + g1$twologlik)
      
      if(traceIt){
        message(sprintf("g=%d, parDiff=%.e, b0=%e, theta0=%e, phi0=%e", 
                        g, parDiff, b0, theta0, phi0))
      }
      
      if(parDiff < 1e-8) break
      
    }
    
    return(list(b=b0, theta=theta0, phi=phi0, logLik=logLik,
         lrt=logLik[length(logLik)] - logLikNULL, 
         lrtASE=2.0*(a1$logLikH1 - a1$logLikH0)))
  }
