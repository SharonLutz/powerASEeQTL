powerASEeQTL <-
function(n,mu=500, n.simu=200,sim="sim1", methods=c("eQTL.linear","eQTL.negBin","eQTL.quasipoisson","ASE", "eQTL.ASE"), folds= seq(1, 2,length.out=3), alpha=0.001, phi=1.1, theta=0.1, maf=0.2, propASE=0.005, legend=TRUE, color=TRUE, title="", subtitle="", titlecolor="black", subtitlecolor="black", titlesize=1, subtitlesize=1, labelsize = 1, labelcolor = "black", linewidth=2,tilt=0, SEED = 1){

  library(VGAM)
  library(MASS)
  
  set.seed(SEED)

  # Check that sim is either sim1 or sim2
  if(!sim%in%c("sim1","sim2")){
    stop("Error: sim must be either sim1 or sim2.")
  }
  
  # Check that n is over 25
  if(n<=25){
    stop("Error: n must be larger.")
  }

possibleMethods <- c("eQTL.linear", "eQTL.negBin","eQTL.quasipoisson","ASE", "eQTL.ASE")

matResultsF<-matrix(0,nrow=length(folds),ncol=length(possibleMethods))
colnames(matResultsF)<-possibleMethods

# Length of methods to be used for color/ line type/ symbol vectors
numMethods <- length(methods)

colorVec <- c("black","red","green3","blue","purple","magenta","mediumorchid4")
colorVec <- c(colorVec[1:numMethods])

# Check that the methods are entered correctly
for(yy in 1:length(methods)){
	if(!toupper(methods[yy])%in%toupper(possibleMethods)){
	  stop(paste("Error: ",methods[yy]," is not in possible methods. See ?powerASEeQTL.",sep=""))
	}
}

### Reorder and rename methods to be same order every time
tempMethods <- rep(NA,length(possibleMethods))

for(jj in 1:length(methods)){
  if(toupper(methods[jj])=="EQTL.LINEAR"){tempMethods[1]<-"eQTL.linear"}
  if(toupper(methods[jj])=="EQTL.NEGBIN"){tempMethods[2]<-"eQTL.negBin"}
  if(toupper(methods[jj])=="EQTL.QUASIPOISSON"){tempMethods[3]<-"eQTL.quasipoisson"}
  if(toupper(methods[jj])=="ASE"){tempMethods[4]<-"ASE"}
  if(toupper(methods[jj])=="EQTL.ASE" ){tempMethods[6]<-"eQTL.ASE"}
}

methods <- tempMethods[!is.na(tempMethods)]
####

# Check that n, mu and n.simu are integers >1
if(!(n>1 & n%%1==0)){stop("Error: n must be an integer > 1.")}
if(!(mu>1 & mu%%1==0)){stop("Error: mu must be an integer > 1.")}
if(!(n.simu>1 & n.simu%%1==0)){stop("Error: n.simu must be an integer > 1.")}

# Check that legend/color = TRUE/FALSE
if(legend!=TRUE & legend!=FALSE){stop("Error: legend must be TRUE or FALSE.")}
if(color!=TRUE & color!=FALSE){stop("Error: color must be TRUE or FALSE.")}

# Check that the min(folds) is greater than or equal to 1
if(min(folds)<1){stop("Error: min(folds) must be >=1.")}

# Check that 0 < alpha < 1 and 0 < maf < 1
if( alpha<=0 | alpha>=1){stop("Error: alpha must be between 0 and 1.")}
if( maf<=0 | maf>=1){stop("Error: maf must be between 0 and 1.")}

# Check that phi > 0 and theta > 0
if(phi<=0){stop("Error: phi must be greater than 0.")}
if(theta<0){stop("Error: theta must be greater than or equal to 0.")}

#Allow for grayscale by creating color vector based on color=TRUE/FALSE
if(color){
  colVec<- colorVec
}else{
  colVec <- gray.colors(numMethods)
}

#Make linetype/symbol vectors
pchv <- c(1:numMethods)
ltyv <- c(1:numMethods)

########
# power simulation
########

for(i in 1:length(folds)){
print(paste("fold change",i,"in",length(folds)))
fold   = folds[i]

matResultsF[i,] = powerRNASEQ(n, mu, fold, phi, theta, n.simu, alpha, maf, methods, sim, propASE)

matResultsR<-matrix(0,nrow=length(folds),ncol=length(methods))
colnames(matResultsR)<-methods

matResultsR[,1]<-matResultsF[,colnames(matResultsF)==methods[1]]
if(length(methods)>1){
	for(gg in 2:length(methods)){
matResultsR[,gg]<-matResultsF[,colnames(matResultsF)==methods[gg]]	
	}}

########
#Plot
########

#matResults<-read.table(file=paste("simulationN",n,"Mu",mu,".txt",sep=""),header=TRUE)

#create pdf of plot
pdf(paste("powerASEeQTLn",n,"Mu",mu,sim,".pdf",sep=""))

plot(-1,-1,xlim=c(min(folds),max(folds)),ylim=c(0,1),main=title,sub=subtitle,col.main=titlecolor,col.sub=subtitlecolor,
     ylab="Power",xlab="Fold Change",col.lab=labelcolor,cex.lab=labelsize,cex.main=titlesize,cex.sub=subtitlesize,las=tilt)

for(mi in 1:length(methods)){

    if(methods[mi]=="eQTL.linear"){lines(folds,matResultsR[,"eQTL.linear"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="eQTL.negBin"){lines(folds,matResultsR[,"eQTL.negBin"],pch = pchv[mi],col=colVec[mi],type="b",lty = ltyv[mi], lwd = linewidth)}

    if(methods[mi]=="eQTL.quasipoisson"){lines(folds,matResultsR[,"eQTL.quasipoisson"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="ASE"&!NA%in%matResultsF[,"ASE"]){lines(folds,matResultsR[,"ASE"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="eQTL.ASE"&!NA%in%matResultsF[,"eQTL.ASE"]){lines(folds,matResultsR[,"eQTL.ASE"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

}

# Create a legend only out of the included methods
leg.text<-rep(0,length(methods))
for(ff in 1:length(methods)){
	if(methods[ff]=="eQTL.linear"){leg.text[ff]<-c("eQTL: Linear Regression")}
	if(methods[ff]=="eQTL.negBin"){leg.text[ff]<-c("eQTL: Negative Binomial")}
	if(methods[ff]=="eQTL.quasipoisson"){leg.text[ff]<-c("eQTL: Quasi Poisson")}
	if(methods[ff]=="ASE"&!NA%in%matResultsF[,"ASE"]){leg.text[ff]<-c("ASE: Beta Binomial")}
  if(methods[ff]=="eQTL.ASE"&!NA%in%matResultsF[,"eQTL.ASE"]){leg.text[ff]<-c("eQTL & ASE: TreCASE")}
}


if(legend==TRUE){legend("bottomright",legend = leg.text,col=colVec,lty=ltyv,pch=pchv)}

}

dev.off()

# Return matResults
matResultsF

}
