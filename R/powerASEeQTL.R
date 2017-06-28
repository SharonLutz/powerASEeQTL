powerASEeQTL <-
function(n,mu=500, n.simu=200, methods=c("linear", "negativeB","poisson","TReC","ASEchi2", "ASEbinom", "TReCASE"), 
         legend=TRUE, color=TRUE, folds= seq(1.5, 1.7, by=0.2), alpha=0.001, phi=1.0, theta=0.1, maf=0.2, title="", subtitle="", 
         titlecolor="black", subtitlecolor="black", titlesize=1, subtitlesize=1, labelsize = 1, labelcolor = "black", linewidth=2, 
         tilt=0, SEED = 1){

  library(VGAM)
  library(MASS)
  
  set.seed(SEED)

possibleMethods <- c("linear", "negativeB", "poisson","TReC","ASEchi2","ASEbinom", "TReCASE")

matResultsF<-matrix(0,nrow=length(folds),ncol=7)
colnames(matResultsF)<-possibleMethods

# Length of methods to be used for color/ line type/ symbol vectors
numMethods <- length(methods)

colorVec <- c("black","red","green3","blue","cyan","magenta","mediumorchid4")
colorVec <- c(colorVec[1:numMethods])

# Check that the methods are entered correctly
for(yy in 1:length(methods)){
	if(!toupper(methods[yy])%in%toupper(possibleMethods)){
	  stop(paste("Error: ",methods[yy]," is not in possible methods. See ?powerASEeQTL.",sep=""))
	}
}

### Reorder and rename methods to be same order every time
tempMethods <- rep(NA,7)

for(jj in 1:length(methods)){
  if(toupper(methods[jj])=="LINEAR"){tempMethods[1]<-"linear"}
  if(toupper(methods[jj])=="NEGATIVEB"){tempMethods[2]<-"negativeB"}
  if(toupper(methods[jj])=="POISSON"){tempMethods[3]<-"poisson"}
  if(toupper(methods[jj])=="TREC"){tempMethods[4]<-"TReC"}
  if(toupper(methods[jj])=="ASECHI2"){tempMethods[5]<-"ASEchi2"}
  if(toupper(methods[jj])=="ASEBINOM"){tempMethods[6]<-"ASEbinom"}
  if(toupper(methods[jj])=="TRECASE"){tempMethods[7]<-"TReCASE"}
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

matResultsF[i,] = powerRNAseq(n, mu, fold, phi, theta, n.simu, alpha, maf)

matResultsR<-matrix(0,nrow=length(folds),ncol=length(methods))
colnames(matResultsR)<-methods

matResultsR[,1]<-matResultsF[,colnames(matResultsF)==methods[1]]
if(length(methods)>1){
	for(gg in 2:length(methods)){
matResultsR[,gg]<-matResultsF[,colnames(matResultsF)==methods[gg]]	
	}}

# outputs the results to the following file in your working directory
write.table(matResultsF, file=paste("simulationN",n,"Mu",mu,".txt",sep=""),quote=F,row.names=F)


########
#Plot
########

#matResults<-read.table(file=paste("simulationN",n,"Mu",mu,".txt",sep=""),header=TRUE)

#create pdf of plot
pdf(paste("powerRNAseqN",n,"Mu",mu,".pdf",sep=""))

plot(-1,-1,xlim=c(min(folds),max(folds)),ylim=c(0,1),main=title,sub=subtitle,col.main=titlecolor,col.sub=subtitlecolor,ylab="Power",xlab="Fold Change",col.lab=labelcolor,cex.lab=labelsize,cex.main=titlesize,cex.sub=subtitlesize,las=tilt)


for(mi in 1:length(methods)){

    if(methods[mi]=="linear"){lines(folds,matResultsR[,"linear"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="negativeB"){lines(folds,matResultsR[,"negativeB"],pch = pchv[mi],col=colVec[mi],type="b",lty = ltyv[mi], lwd = linewidth)}

    if(methods[mi]=="poisson"){lines(folds,matResultsR[,"poisson"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="TReC"){lines(folds,matResultsR[,"TReC"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="ASEchi2"){lines(folds,matResultsR[,"ASEchi2"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="ASEbinom"){lines(folds,matResultsR[,"ASEbinom"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

    if(methods[mi]=="TReCASE"){lines(folds,matResultsR[,"TReCASE"],pch = pchv[mi],col = colVec[mi],type="b",lty = ltyv[mi],lwd = linewidth)}

}

# Create a legend only out of the included methods
leg.text<-rep(0,length(methods))
for(ff in 1:length(methods)){
	if(methods[ff]=="linear"){leg.text[ff]<-c("eQTL: Linear Regression")}
	if(methods[ff]=="negativeB"){leg.text[ff]<-c("eQTL: Negative binomial")}
	if(methods[ff]=="poisson"){leg.text[ff]<-c("eQTL: Poisson Regression")}
  if(methods[ff]=="TReC"){leg.text[ff]<-c("eQTL:TreC")}
  if(methods[ff]=="ASEchi2"){leg.text[ff]<-c("ASE: chi-square")}
  if(methods[ff]=="ASEbinom"){leg.text[ff]<-c("ASE: binomial test")}
  if(methods[ff]=="TReCASE"){leg.text[ff]<-c("eQTL & ASE: TreCASE")}
}

if(legend==TRUE){legend("topleft",legend = leg.text,fill=colVec,col=colVec,lty=ltyv,pch=pchv)}

}

dev.off()

# Return matResults
matResultsF

}
