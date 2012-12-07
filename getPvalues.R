########################################################################################
# Functions to scrape P-values from Pubmed abstracts
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
#          read on windows machines. 
#
########################################################################################



library(RCurl)
library(XML)
library(tm)


# A function to get all abstracts and pubmed ids for papers from the journal "journaltitle" in the year "year"
# by scraping the Pubmed API.

getAbstractsPmids = function(journaltitle,year){
# esearch
url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
q = paste("db=pubmed&term=",gsub(" ","+",journaltitle),"[ta]+AND+",year,"[dp]&usehistory=y",sep="")
esearch <- xmlTreeParse(getURL(paste(url, q, sep="")), useInternal = T)
webenv  <- xmlValue(getNodeSet(esearch, "//WebEnv")[[1]])
key     <- xmlValue(getNodeSet(esearch, "//QueryKey")[[1]])

# efetch
url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
q   <- "db=pubmed&retmode=xml&rettype=abstract"
efetch <- xmlTreeParse(getURL(paste(url, q, "&WebEnv=", webenv, "&query_key=", key, sep="")), useInternal = T)
r = xmlRoot(efetch)

n = xmlSize(r)
abstracts = pmid = titles = rep(NA,n)
for(i in 1:n){abstracts[i] =  xmlValue(r[[i]][[1]][["Article"]][["Abstract"]]); pmid[i] = xmlValue(r[[i]][[1]][["PMID"]]); titles[i] = xmlValue(r[[i]][[1]][["Article"]][["ArticleTitle"]]) }
return(list(abstracts=abstracts,pmid=pmid,titles=titles))
}


#A function to remove trailing zeros from the P-value strings
removeTrailing = function(string){
  while(length(grep("[0-9]",strsplit(string,"")[[1]][nchar(string)])) == 0){
    string = substr(string,1,(nchar(string)-1))
  }
  return(string)
}


# A function to convert the scientific notation used by journals into
# numeric values that can be analyzed. 
convertScientific = function(string){
  if(length(grep("[[:punct:]][[:space:]]",string))>0){
    string = strsplit(string,"[[:punct:]][[:space:]]")[[1]][1]
    string = removeTrailing(string)
  }
  if(length(grep("[×x]",string))>0){
    string = gsub("[:space:]","",string)
    tmp1 = as.numeric(strsplit(string,"[×x]")[[1]][1])
    tmp2 = as.numeric(strsplit(strsplit(string,"[×x]")[[1]][2],"[−-]")[[1]][2])
    return(tmp1*10^(-tmp2))
  }else{
    return(as.numeric(string))
  }
}


# A function to scrape the P-values from a vector of abstracts with corresponding
# pubmed ids

getPvalues = function(abstract,pmid){
  pvalues = numeric(0)
  trunc = numeric(0)
  ids = numeric(0)

  # Get the truncated p-values
  ind = grep("[Pp][[:space:]]?[<≤]",abstract)
  for(i in 1:length(ind)){
    tmp = strsplit(abstract[ind[i]],"[[:space:](][Pp][[:space:]]?[<≤]")[[1]]
    n = length(tmp)
    for(j in 1:n){
      if(length(grep("[.0123456789]",substr(tmp[j],1,2))) > 0){
        if(length(grep("[A-Z]",substr(tmp[j],1,1)))>0){next;}
        tmp2 = strsplit(tmp[j],"[^[:punct:][:digit:]x[:space:]]")[[1]][1]
        tmp2 = removeTrailing(tmp2)
        tmp2 = gsub(" ","",tmp2)
        tmp2 = convertScientific(tmp2)
        pvalues = c(pvalues,as.numeric(tmp2))
        trunc = c(trunc,1)
        ids = c(ids,pmid[ind[i]])
      }
    }
  }

   # Get the truncated p-values
  ind = grep("[Pp][[:space:]]?=",abstract)
  for(i in 1:length(ind)){
    tmp = strsplit(abstract[ind[i]],"[[:space:](][Pp][[:space:]]?=")[[1]]
    n = length(tmp)
    for(j in 1:n){
      if(length(grep("[.0123456789]",substr(tmp[j],1,2))) > 0){
        if(length(grep("[A-Z]",substr(tmp[j],1,1)))>0){next;}
        tmp2 = strsplit(tmp[j],"[^[:punct:][:space:][:digit:]x[:space:]]")[[1]][1]
        tmp2 = removeTrailing(tmp2)
        tmp2 = gsub(" ","",tmp2)
        tmp2 = convertScientific(tmp2)
        pvalues = c(pvalues,as.numeric(tmp2))
        trunc = c(trunc,0)
        ids = c(ids,pmid[ind[i]])
      }
    }
  }
  return(list(pvalues=pvalues,ids=ids,trunc=trunc))
}


## Run the above functions over a specified set of journals and years
## to obtain P-values


journals = c("JAMA","New England Journal of Medicine","BMJ","American Journal of Epidemiology","Lancet")
years = 2000:2010

pvalueData = matrix(NA,nrow=1,ncol=6)
colnames(pvalueData) = c("pvalue","pvalueTruncated","pubmedID","year","abstract","title")

npapers = matrix(NA,nrow=length(journals),ncol=length(years))

for(i in 1:length(journals)){
  for(j in 1:length(years)){
    cat(journals[i]); cat(" "); cat(years[j]); cat(" "); 
    tmpData = getAbstractsPmids(journals[i],years[j])
    while(length(tmpData$abstracts) ==1 & is.na(tmpData$abstracts[1])){tmpData = getAbstractsPmids(journals[i],years[j])}
    cat("Downloaded"); cat(" ");
    npapers[i,j] = length(tmpData$abstracts)
    tmpOut = getPvalues(tmpData$abstracts,tmpData$pmid)
    
    nPvalues = length(tmpOut$pvalues)
    aa = match(tmpOut$ids,tmpData$pmid)
    
    tmpMatrix = cbind(tmpOut$pvalues,tmpOut$trunc,as.numeric(tmpOut$ids),rep(years[j],nPvalues),tmpData$abstracts[aa],tmpData$titles[aa])
    rownames(tmpMatrix) = rep(journals[i],nPvalues)
    pvalueData = rbind(pvalueData,tmpMatrix)
    cat("Done\n")
  }
}


## These come from the paper Pubmed ID = 20876667 and are reports from
## a GWAS. Update some P-values that are incorrectly scraped because of
## variation in scientific notation.

pvalueData[pvalueData[,1] == 10,1] = 10e-7
pvalueData = pvalueData[-1,]

## These come from the Lancet, with missed periods in the P-values

pvalueData[13392,1] = 0.003
pvalueData[13413,1] = 0.014
pvalueData[13414,1] = 0.004

## This one comes from the Lancet too, where an incorrect P-value is grabbed.

pvalueData = pvalueData[-14674,]


# Remove rows with P-values that are 
pvalueData = pvalueData[!is.na(pvalueData[,1]),]

# This one for some reason replaced a period with "small middle dot"

pvalueData[which(pvalueData[,3] == 11943262),1] = 1e-4

# This one has a <<< in the P-value definition
pvalueData = pvalueData[-which(pvalueData[,3]=="16421237"),] 


save(pvalueData,npapers,file="pvalueData.rda")



