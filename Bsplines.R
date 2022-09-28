#!/usr/bin/env Rscript
library(dplyr)
library(splines)
library(tidyverse)

normalize <- function(x){
  return((x-min(x)) / (max(x)-min(x)))
}

args = commandArgs(trailingOnly=TRUE)
print(args[1])
samplename <- str_split(args[1], "avg")[[1]][1]
data <- read.table(args[1], sep='\t')
X<-data%>%select(colnames(data))%>%as.matrix()
X <- normalize(X)
#print(X)
x=seq(0,1,length.out=ncol(X))
B = bs(x, df=7, degree = 3)[,1:7]
Bcoef = matrix(0,dim(X)[1],7)

for(i in 1:dim(X)[1])
{
Bcoef[i,] = solve(t(B)%*%B)%*%t(B)%*%as.matrix(X[i,])
}
Bcoef <- data.frame(Bcoef, genes=rownames(data))
colnames(Bcoef)[colnames(Bcoef) != "genes"] <-  paste(samplename, colnames(Bcoef)[colnames(Bcoef) != "genes"], sep="_")
write.table(Bcoef,paste0(args[1],'bsplines.tsv'), sep='\t')
