
powerASEeQTL <-
function(n,mu=500, n.simu=200, methods=c("Linear Model", "Negative Binomial", "TReC", "TReCASE", "ASEchi2","ASEbinom","Poisson"), legend=TRUE, color=TRUE, folds= seq(1.5, 1.7, by=0.2), alpha=0.001, phi=1.0, theta=0.1, maf=0.2, numcol=1, title="", subtitle="", titlecolor="black", subtitlecolor="black", titlesize=1, subtitlesize=1, legendbox="o", labelsize = 1, labelcolor = "black", linewidth=2, tilt=0)
{
  
  library(VGAM)
  library(MASS)

possibleMethods <- c("Linear Model", "Negative Binomial", "TReC", "TReCASE", "ASEchi2", "ASEbinom", "Poisson")

numMethods <- length(possibleMethods)


#CHECKS
#Check that the length of folds is > 1
stopifnot(length(folds)>1)
#Check that methods is the right length
stopifnot(length(methods)>=1 & length(methods)<=numMethods)
#Check that the methods are entered correctly
for(i in 1:length(methods)){
	stopifnot(methods[i] %in% possibleMethods)
}
 
#Allow for grayscale by creating color vector based on color=TRUE/FALSE
if(color){
  colVec<- c(1:numMethods)
}else{
  colVec <- gray.colors(numMethods)
}

#Make linetype/symbol vectors
pchv <- c(1:numMethods)
ltyv <- c(1:numMethods)

########
# power simulation
########

matResults<-matrix(0,nrow=length(folds),ncol=numMethods)


colnames(matResults)<-c("Linear","glmNB","TReC","TReCASE","ASEchi2","ASEbinom","Poisson")

for(i in 1:length(folds)){
print(paste("folds",i,"in",length(folds)))

fold   = folds[i]

matResults[i,] = powerRNAseq(n, mu, fold, phi, theta, n.simu, alpha, maf)

#outputs the results to the following file in your working directory
write.table(matResults, file=paste("simulationN",n,"Mu",mu,".txt",sep=""),quote=F,row.names=F)


########
#Plot
########

matResults<-read.table(file=paste("simulationN",n,"Mu",mu,".txt",sep=""),header=TRUE)


#create pdf of plot
pdf(paste("powerRNAseqN",n,"Mu",mu,".pdf",sep=""))


plot(-1,-1,xlim=c(min(folds),max(folds)),ylim=c(0,1),main=title,sub=subtitle,col.main=titlecolor,col.sub=subtitlecolor,ylab="Power",xlab="Fold Change",col.lab=labelcolor,cex.lab=labelsize,cex.main=titlesize,cex.sub=subtitlesize,las=tilt)


for(mi in 1:length(methods)){
    if(methods[mi]=="Linear Model"){lines(folds,matResults[,1],pchv[mi],colVec[mi],type="b",ltyv[mi],lwd=linewidth)}
    if(methods[mi]=="Negative Binomial"){lines(folds,matResults[,2],pchv[mi],colVec[mi],type="b",ltyv[mi], lwd=linewidth)}
    if(methods[mi]=="TReC"){lines(folds,matResults[,3],pchv[mi],colVec[mi],type="b",ltyv[mi],lwd=linewidth)}
    if(methods[mi]=="TReCASE"){lines(folds,matResults[,4],pchv[mi],colVec[mi],type="b",ltyv[mi],lwd=linewidth)}
    if(methods[mi]=="ASEchi2"){lines(folds,matResults[,5],pchv[mi],colVec[mi],type="b",ltyv[mi],lwd=linewidth)}
    if(methods[mi]=="ASEbinom"){lines(folds,matResults[,6],pchv[mi],colVec[mi],type="b",ltyv[mi],lwd=linewidth)}
    if(methods[mi]=="Poisson"){lines(folds,matResults[,7],pchv[mi],colVec[mi],type="b",ltyv[mi],lwd=linewidth)}
}

if(legend){legend(min(folds),1,methods,fill=colVec,col=colVec,lty=ltyv,pch=pchv,ncol=numcol,bty=legendbox)}

}

dev.off()


}