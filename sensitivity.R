
##
## Load relevant libraries
##

library(RColorBrewer)
source("calculateSwfdr.R")
source("journalAnalysisHelp.R")
library(stats4)


##
## Functions for creating rounded and truncated p-values
##

roundTrunc = function(pp,probTrunc,probRound){
  tt = rbinom(length(pp),prob=probTrunc,size=1)
  rr = rbinom(length(pp),prob=probRound,size=1)
  rr[tt==1 & rr==1] <- 0
  pp[rr==1] = round(pp[rr==1],2)
  pp[tt==1] = truncPvals(pp[tt==1])
  return(list(pp=pp,tt=tt,rr=rr))
}

truncPvals = function(pp){
  pp[pp < 0.001] <- 0.001
  pp[pp > 0.001 & pp < 0.01] <-  0.01
  pp[pp > 0.01 & pp < 0.02] <- 0.02
  pp[pp > 0.02 & pp < 0.03] <- 0.03
  pp[pp > 0.03 & pp < 0.04] <- 0.04
  pp[pp > 0.04 & pp <= 0.05] <- 0.05
  return(pp)
}



###
### Simulate cases where the alternative distribution stays
### fixed, but the null p-value distribution is more and more
### anti-conservative (right skewed)
###

set.seed(8383)
ntotal = 1000
fdr = seq(0.1,0.9,by=0.1)
bvalue = c(1,5,10,25,50,75,100)
fdrhat = matrix(0,ncol=length(fdr),nrow=length(bvalue))

for(j in 1:length(bvalue)){
  for(i in 1:9){
    # Create p-vals
    nalt = ntotal*(1-fdr[i])
    nnull = ntotal*fdr[i]
    altpvals = rbeta(1e5,1,100)
    altpvals = altpvals[altpvals < 0.05][1:nalt]
    nullpvals = rbeta(1e5,1,bvalue[j])
    nullpvals = nullpvals[nullpvals < 0.05][1:nnull]
    pvals = c(altpvals,nullpvals)
    tmp = roundTrunc(pvals,0.1,0.1)
    fdrhat[j,i] = calculateSwfdr(tmp$pp,tmp$tt,tmp$rr,numEmIterations=100)$pi0
    cat(i)
  }
  print(j)
}


pdf(file="sensitivity.pdf",height=4,width=8)
mypar <- function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
 par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
 par(mfrow=c(a,b),...)
 palette(brewer.pal(brewer.n,brewer.name))
}



mypar(mfrow=c(1,2))
plot(1:10,1:10,type="n",xlab="SWFDR",ylab="Estimated SWFDR",xlim=c(0,1),ylim=c(0,1))
for(i in 1:length(bvalue)){
lines(fdr,fdrhat[i,],xlim=c(0,1),ylim=c(0,1),xlab="True swfdr",ylab="Estimated swfdr",pch=22,bg=i,col="black",cex.axis=1.5,cex.lab=1.5,cex=2,type="b")
}
abline(c(0,1),lwd=3,col="grey ")


xvals = seq(0,0.05,length=100)
plot(xvals,dbeta(xvals,1,100)/pbeta(0.05,1,100),col="black",type="l",lwd=2,xlab="P-value",ylab="Density")
for(i in 1:length(bvalue)){
  lines(xvals,dbeta(xvals,1,bvalue[i])/pbeta(0.05,1,bvalue[i]),col=i,lwd=3,lty=2)
}
dev.off()



