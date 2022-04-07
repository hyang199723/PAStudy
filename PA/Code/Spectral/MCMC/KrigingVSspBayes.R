# Compare our Kriging algorithm result with spBayes
rm(list = ls())
library(fields) 
library(geoR)
library(truncnorm)
library(tidyverse)
library(mvtnorm)
library(gridExtra)
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Spectral/MCMC/")
source('ExtraFunctions.R')
source('LMC_function.R')

#simulate data
source('simAllTS.R') # load your data here 

Y11 = Y1[train.sites1,]
s11=s1[train.sites1,]
sp1=s1[which.test1,]


# train set of PA (all)
Y22 = Y2
s22=s2
sp2=NULL # no test set

n1       <- nrow(Y11)
n2       <- nrow(Y22)
nt       <- ncol(Y11)
m1       <- is.na(Y11)
m2       <- is.na(Y22)
d        <- as.matrix(dist(rbind(s11,s22)))
dv2 = as.matrix(dist(s22))

np1 = nrow(sp1)
np2  = nrow(sp2)
if (is.null(sp2))
{np2  = 0}
dp = as.matrix(dist(rbind(sp1,sp2)))
all.d=as.matrix(dist(rbind(s1,s2,sp1,sp2))) 
Z1p = matrix(0,np1,nt)


# Kriging 
Mp=exp_corr(all.d,range=rangeu)
Mp00=Mp[1:(n1+n2),1:(n1+n2)]
Mp11=Mp[(n1+n2+1):(n1+n2+np1+np2),(n1+n2+1):(n1+n2+np1+np2)]
Mp10=Mp[(n1+n2+1):(n1+n2+np1+np2),1:(n1+n2)]
Mp01=t(Mp10)

E00=eigen(Mp00)
E00.G=E00$vectors
E00.D=E00$values

Mp00.inv=E00.G%*%diag(1/E00.D)%*%t(E00.G)

AA=Mp10%*%Mp00.inv
a=Mp10%*%Mp00.inv%*%t(Mp10)
a=round(a,digits=7) # to avoid numercial underflow
B=Mp11-a

### sample U npreds times 
for (j in 1:npreds)
{
  Uls = matrix(0,n1+n2,nt)
  
  for (r in 1:nt)
  {
    Uls[,r]=t(chol(Mp00))%*%rnorm(n1+n2,0,sqrt(sigmau[r]))
  }
  
  U1p = matrix(0,np1,nt)
  for (r in 1:nt)
  {
    Au=AA%*%Uls[,r] 
    sigmaB=sigmau[r]*B
    Ul.pred=rmvnorm(1,mean=Au,sigma=sigmaB)
    U1p[,r]=Ul.pred[(1:np1)]
  }
  
  for(nn in 1:np1){Z1p[nn,] <- fft_real(U1p[nn,],inverse=TRUE)}
  Y1.pred[which.test1,,j] <- beta.1+Z1p+rnorm(n=nt,sd=sqrt(tau1)) 
}