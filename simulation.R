###
### Simulate 100 different 'journals' where the p-values
### that are reported are all the p-values less than 0.05
###


### Source our functions
###

source("calculateSwfdr.R")
source("journalAnalysisHelp.R")
library(stats4)


##
## Perform the simulation
##

set.seed(23452134)

n=100
pi0 = pi0hat = rep(NA,n)
for(i in 1:n){

palt = rbeta(1e3,1,100)
palt = palt[palt < 0.05]
x = runif(1,0,5)
pnull = runif(floor(length(palt)*x),0,0.05)

pp = c(pnull,palt)



tt = rbinom(length(pp),size=1,prob=0.2)
rr = rbinom(length(pp),size=1,prob=0.2)
rr[tt==1] = 0

pp[rr > 0] = round(pp[rr > 0],2)


pi0[i] = length(pnull)/length(pp)
pi0hat[i] = calculateSwfdr(pp,tt,rr,alpha=1,beta=150,numEmIterations=10)$pi0
print(i)
}

### Plot our estimate of the swfdr versus the truth. 
###

pdf(file="all-significant.pdf")
plot(pi0,pi0hat,xlim=c(0,1),ylim=c(0,1),xlab="True swfdr",ylab="Estimated swfdr",pch=22,bg="blue",col="black",cex.axis=1.5,cex.lab=1.5,cex=2)
abline(c(0,1),lwd=3,col="grey ")
dev.off()






###
### Simulate from the case where only the minimum
### p-value is reported if it less than 0.05
###

n=100
pi0 = pi0hat = rep(NA,n)
for(i in 1:n){

palt = rep(NA,1000)
for(j in 1:1000){
 ptmp = rbeta(20,1,100)
 pval = min(ptmp)
 if(pval < 0.05){palt[j] = pval}
}
palt = palt[!is.na(palt)]


x = runif(1,0,5)
pnull = rep(NA,floor(length(palt)*x))

for(j in 1:length(pnull)){
  ptmp = rbeta(20,1,1)
  pval = min(ptmp)
  if(pval < 0.05){pnull[j] = pval}
}
pnull = pnull[!is.na(pnull)]

pp = c(pnull,palt)

tt = rbinom(length(pp),size=1,prob=0.2)
rr = rbinom(length(pp),size=1,prob=0.2)
rr[tt==1] = 0

pp[rr > 0] = round(pp[rr > 0],2)


pi0[i] = length(pnull)/length(pp)
pi0hat[i] = calculateSwfdr(pp,tt,rr,alpha=1,beta=150,numEmIterations=10)$pi0
print(i)
}

pdf(file="only-min.pdf")
plot(pi0,pi0hat,xlim=c(0,1),ylim=c(0,1),xlab="True swfdr",ylab="Estimated swfdr",pch=22,bg="blue",col="black",cex.axis=1.5,cex.lab=1.5,cex=2)
abline(c(0,1),lwd=3,col="grey")
dev.off()



###
### Simulate from the case where p-values are always
### rounded down when they are rounded. 
###

n=100
pi0 = pi0hat = rep(NA,n)
for(i in 1:n){

palt = rep(NA,1000)
for(j in 1:1000){
 ptmp = rbeta(20,1,100)
 pval = min(ptmp)
 if(pval < 0.05){palt[j] = pval}
}
palt = palt[!is.na(palt)]


x = runif(1,0,5)
pnull = rep(NA,floor(length(palt)*x))

for(j in 1:length(pnull)){
  ptmp = rbeta(20,1,1)
  pval = min(ptmp)
  if(pval < 0.05){pnull[j] = pval}
}
pnull = pnull[!is.na(pnull)]

pp = c(pnull,palt)

tt = rbinom(length(pp),size=1,prob=0.2)
rr = rbinom(length(pp),size=1,prob=0.2)
rr[tt==1] = 0

pp[rr > 0] = floor(pp[rr > 0]*100)/100


pi0[i] = length(pnull)/length(pp)
pi0hat[i] = calculateSwfdr(pp,tt,rr,alpha=1,beta=150,numEmIterations=10)$pi0
print(i)
}

pdf(file="floored-rounding.pdf")
plot(pi0,pi0hat,xlim=c(0,1),ylim=c(0,1),xlab="True swfdr",ylab="Estimated swfdr",pch=22,bg="blue",col="black",cex.axis=1.5,cex.lab=1.5,cex=2)
abline(c(0,1),lwd=3,col="grey")
dev.off()
