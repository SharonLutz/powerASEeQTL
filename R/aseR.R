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
  
  list(parH0=theta0,  logLikH0=op0$value, 
       parH1=op1$par, logLikH1=op1$value)
}
