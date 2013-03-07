########################################################################################
# Analysis of P-value data to estimate journal level false positive rates
# Date: 7-1-12
# Copyright (C) 2011 Jeffrey T. Leek (http://www.biostat.jhsph.edu/~jleek/contact.html) and Leah R. Jager
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details, see <http://www.gnu.org/licenses/>.
#
#
#    Note: These functions were written on a Mac and may have difficulties when
#          read on Windows machines. They depend on the files "calculateSwfdr.R" and
#          "journalAnalysisHelp.R" and the pvalue data calculated with "getPvalues d.R"
#           all code is available from: https://github.com/jtleek/swfdr
#          It also depends on the R libraries: stat4, genefilter, lme4
#
########################################################################################



###
### Section 1 - Calculate quantities of interest
###

### Load the data and functions

load("pvalueData.rda")
cNames = colnames(pvalueData)[1:4]
rNames = rownames(pvalueData)
pvalueData = matrix(as.numeric(pvalueData[,1:4]),ncol=4)
colnames(pvalueData) = cNames
rownames(pvalueData) = rNames

source("calculateSwfdr.R")
source("journalAnalysisHelp.R")
library(stats4)
library(genefilter)
library(lme4)

# Initialize the variables
percent05 = pi0JournalYear = matrix(NA,nrow=5,ncol=11)
pi0Journal = rep(NA,5)

# Get the journals and years
journals = c("Lancet", "JAMA","New England Journal of Medicine","BMJ","American Journal of Epidemiology")
journalAbb =c("Lancet","JAMA","NEJM","BMJ","AJE")
years <- 2000:2010

# Get the percent of P-values less than 0.05 by journal/year

for(i in 1:length(journals)){
  for(j in 1:length(years)){
    percent05[i,j] = mean(as.numeric(pvalueData[pvalueData[,4]==years[j] & rownames(pvalueData) == journals[i],1]) < 0.05,na.rm=T)
  }
}



# Estimated the swfdr by journal/year

for(i in 1:length(journals)){
  for(j in 1:length(years)){
    tmp = pvalueData[pvalueData[,4]==years[j] & rownames(pvalueData) == journals[i] & pvalueData[,1] < 0.05,]
    tt = tmp[,2]
    rr = rep(0,length(tt))
    rr[tt == 0] = (tmp[tt==0,1] == round(tmp[tt==0,1],2))
    pp = tmp[,1]
    out = calculateSwfdr(pValues = pp, truncated = tt, rounded = rr, numEmIterations=100)
    pi0JournalYear[i,j] = out$pi0
    cat(j)
  }
  cat(i)
  cat("\n")
}

zlist = plist = rlist = tlist = vector(length(journals),mode="list")

# Estimate the swfdr by journal
for(i in 1:length(journals)){
  
  tmp = pvalueData[rownames(pvalueData) == journals[i] & as.numeric(pvalueData[,1]) < 0.05,]
  tt = as.numeric(tmp[,2])
  rr = rep(0,length(tt))
  rr[tt == 0] = (as.numeric(tmp[tt==0,1]) == round(as.numeric(tmp[tt==0,1]),2))
  pp = tmp[,1]
  out = calculateSwfdr(pValues = pp, truncated = tt, rounded = rr, numEmIterations=100)
  zlist[[i]] = out$z
  plist[[i]] = pp
  rlist[[i]] = rr
  tlist[[i]] = tt
  pi0Journal[i] = out$pi0
  cat(i)
}

journalAbb
pi0Journal



pi0JournalBoot = matrix(NA,nrow=5,ncol=100)
# Get the boostrapped statistics by journal
set.seed(839283)
for(j in 1:100){
  for(i in 1:length(journals)){
    tmp = pvalueData[rownames(pvalueData) == journals[i] & pvalueData[,1] < 0.05,]
    tmp = tmp[sample(1:dim(tmp)[1],size=dim(tmp)[1],replace=TRUE),]
    tt = tmp[,2]
    rr = rep(0,length(tt))
    rr[tt == 0] = (tmp[tt==0,1] == round(tmp[tt==0,1],2))
    pp = tmp[,1]
    out = calculateSwfdr(pValues = pp, truncated = tt, rounded = rr, numEmIterations=100)
    pi0JournalBoot[i,j] = out$pi0
    cat(i)
  }
}

apply(pi0JournalBoot,1,sd)

save(pi0JournalBoot,file="pi0JournalBoot.rda")

# Estimate the overall swfdr

  
  tmp = pvalueData[as.numeric(pvalueData[,1]) < 0.05,]
  tt = as.numeric(tmp[,2])
  rr = rep(0,length(tt))
  rr[tt == 0] = (as.numeric(tmp[tt==0,1]) == round(as.numeric(tmp[tt==0,1]),2))
  pp = as.numeric(tmp[,1])
  out = calculateSwfdr(pValues = pp, truncated = tt, rounded = rr, numEmIterations=100)
  pi0Overall = out$pi0


pi0Overall

# Get bootstrapped statistics for standard errors

set.seed(839283)
pi0OverallBoot = rep(NA,100)
for(i in 1:100){
  tmp = pvalueData[pvalueData[,1] < 0.05,]
  tmp = tmp[sample(1:dim(tmp)[1],size=dim(tmp)[1],replace=TRUE),]
  tt = tmp[,2]
  rr = rep(0,length(tt))
  rr[tt == 0] = (tmp[tt==0,1] == round(tmp[tt==0,1],2))
  pp = tmp[,1]
  out = calculateSwfdr(pValues = pp, truncated = tt, rounded = rr, numEmIterations=100)
  pi0OverallBoot[i] = out$pi0
  cat(i)
}

sd(pi0OverallBoot)

save(pi0OverallBoot,file="pi0OverallBoot.rda")


numSubmissions = matrix(NA,nrow=length(journals),ncol=11)
rownames(numSubmissions) = journals
numSubmissions[2,] = c(4366, 4189, 4615, 5064, 5184, 5744, 5352, 5551, 5525, 5946, 5878)
numSubmissions[3,] = c(3566,3758,4238,4486,4870,5066,5037,5412,5909,5642,5284)
numSubmissions[5,] = c(649,609,652,705,817,970,947,922,1010,985,1030)
numSubmissions[4,] = c(NA,NA,NA,3477,3738,4051,3398,3207,3361,3792,3486)
numSubmissions[1,] = c(NA,NA,NA,NA,NA,NA,5254,4770,4720,4609,4612)



###
### Section 2 - Perform Analysis 
###

# Get the number of unique Pubmed IDs

length(unique(pvalueData[,3]))

# Get the percent less than 0.05 by journal


percent05Journal <- rep(NA,length(journals))
for(i in 1:length(journals)){
   percent05Journal[i] = mean(as.numeric(pvalueData[rownames(pvalueData) == journals[i],1]) < 0.05,na.rm=T)
}

journalAbb
percent05Journal


# Calculate the rate of increase in false positives by year and number of submissions
journalInd = as.factor(rep(1:5,each=11))
yearsVals = rep(2000:2010,5)
numSubVec = as.vector(t(numSubmissions))
pi0JournalYearVec = as.vector(t(pi0JournalYear))

glmm1 = glmer(pi0JournalYearVec ~ yearsVals + (1|journalInd),family="gaussian")
p.values.lmer(glmm1)

glmm2 = glmer(pi0JournalYearVec ~ numSubVec + (1|journalInd),family="gaussian")
p.values.lmer(glmm2)


# Calculate the rate of increase in false positives by year and number of submissions, journal as a fixed effect

journalInd = as.factor(rep(1:5,each=11))
yearsVals = rep(2000:2010,5)
numSubVec = as.vector(t(numSubmissions))
pi0JournalYearVec = as.vector(t(pi0JournalYear))

lm1 <- lm(pi0JournalYearVec ~as.factor(journalInd) + yearsVals)
summary(lm1)

lm2 <- lm(pi0JournalYearVec ~ as.factor(journalInd) + numSubVec)
summary(lm2)




# Make the figure showing the increasing number of submissions by journals


pdf(file="figure2.pdf",height=8,width=8)
cols = c("blue","red","orange","green","purple","darkgrey")
cx = 3
plot(1:10,1:10,type="n",xlim=c(2000,2010), ylim=c(0,6500),xlab="",ylab="Number of Submissions",cex.lab=1.5,cex.axis=1.5,xaxt="n")
for(i in 1:5){
  points(2000:2010,numSubmissions[i,],bg=cols[i],col="black",cex=cx,type="b",pch=22,lwd=2)

}
legend(2001,3000,legend=journalAbb,fill=cols,cex=1.5)
axis(1,at=years,labels=years,las=2,cex.axis=2)

dev.off()


# Get the mean number of P-values less than 0.05 across all journals

percent05journals = rep(NA,length(journals))
for(i in 1:length(journals)){
    percent05journals[i] = mean(as.numeric(pvalueData[rownames(pvalueData) == journals[i],1]) < 0.05,na.rm=T)
}


## Get the median number of papers per journal

median(table(pvalueData[,3]))
mad(table(pvalueData[,3]))



## Select 10 papers for validation

set.seed(43888431)
x = sample(unique(pvalueData[,3]),replace=F,size=10)
x



### Make Figure 3

cols = c("blue","red","orange","green","purple","darkgrey")


years3 = c(2000,2005,2010)

pdf("figure3.pdf",height=5*8,width=3*8)
par(mfcol=c(5,3),mar=3*c(2,3,1,1))
for(i in 1:3){
  for(j in 1:length(journals)){
    hist(pvalueData[pvalueData[,4]== years3[i] & rownames(pvalueData)==journals[j] & pvalueData[,1] < 0.05,1],breaks=50,freq=T,xlab="",ylab="",xlim=c(0,0.05),main="",col=cols[j],ylim=c(0,210),border="black",cex.axis=2)
    text(0.025,180,paste(journalAbb[j],years3[i]),cex=5)
    mtext("p-value",side=1,line=3,cex=2)
    mtext("Density",side=2,line=3,cex=2)
  }
}
dev.off()


### Make Figure 4


pdf(file="figure4.pdf",height=8,width=8)
cols = c("blue","red","orange","green","purple","darkgrey")
cx = 3
plot(1:10,1:10,type="n",xlim=c(2000,2010), ylim=c(0,1),xlab="",ylab="False Discovery Rate",cex.lab=1.5,cex.axis=1.5,xaxt="n")
for(i in 1:5){
  points(2000:2010,pi0JournalYear[i,],bg=cols[i],col="black",cex=cx,type="b",pch=22,lwd=2)

}
legend(2001,0.9,legend=journalAbb,fill=cols,cex=1.5)
axis(1,at=years,labels=years,las=2,cex.axis=2)

dev.off()






