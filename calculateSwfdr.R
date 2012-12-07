########################################################################################
# Function for calculating the science-wise FDR 
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
#          read on Windows machines. They depend on the functions in "journalAnalysisHelp.R" 
#           all code is available from: https://github.com/jtleek/swfdr
#          It also depends on the R libraries: stat4, genefilter, lme4
#
########################################################################################

calculateSwfdr = function(pValues,truncated,rounded,pi0 = 0.5,alpha=1,beta=50,numEmIterations=100){
  pp = pValues
  tt = truncated
  rr = rounded
  
  
  ll = function(a,b){
    tmp1 = rep(0,length(pp))
    tmp1[tt==0 & rr==0] = log(dbeta(pp[tt==0 & rr==0],a,b)/pbeta(0.05,a,b))
    tmp1[tt > 0 & rr==0] = log(pbeta(pp[tt > 0 & rr==0],a,b)/pbeta(0.05,a,b))
    tmp1 = -sum((1-z)*tmp1)
    
    probvec = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),a,b) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),a,b))/pbeta(0.05,a,b)
    probvec = rev(probvec)
    tmp2 = sum(-n1*log(probvec))
    
  
    return(tmp1 + tmp2)
  }
  
  n = table(cut(pp[rr > 0],c(-0.01,0.005,0.015,0.025,0.035,0.045,0.051)))
  
  
  for(i in 1:numEmIterations){
    
    ## E-step
    
     probvec1 = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),alpha,beta) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),alpha,beta))/pbeta(0.05,alpha,beta)
     probvec1 = rev(probvec1)
     probvec0 = c(0.005,0.01,0.01,0.01,0.01,0.005)*20
     
     pij0 = pi0*probvec0/(probvec0*pi0 + probvec1*(1-pi0))
     n0 = n*pij0
     n1 = n - n0
     
     z = rep(0,length(pp))
     z[tt == 0 & rr ==0] <- pi0*20/(pi0*20 + (1-pi0)*dbeta(pp[tt==0 & rr == 0],alpha,beta)/pbeta(0.05,alpha,beta))
     z[tt > 0 & rr ==0] <- pi0*20*pp[tt > 0 & rr ==0]/(pi0*20*pp[tt > 0 & rr==0] + (1-pi0)*pbeta(pp[tt > 0 & rr==0],alpha,beta)/pbeta(0.05,alpha,beta))
     
     ## M-step
     
     pi0 = (sum(n0) + sum(z))/(sum(n) + sum(rr == 0))
     tmp = mle(ll,start=list(a=0.05,b=100),lower=c(0.001,1),upper=c(1,500),method="L-BFGS-B")
     alpha =coef(tmp)[1]
     beta = coef(tmp)[2]
   }
  return(list(pi0 = pi0, alpha=alpha, beta = beta, z=z,n0=n0,n=n))
}





