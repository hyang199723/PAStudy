# Model Simulation 

library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)

# set some parameters
n=c(25,15) # number of locations


tau1=1 # error variance1
tau2=1.5  # error variance2

# correlation parameters
sigma1=1 
sigma2=4
range1=exp(0.1)
range2=exp(0.3)

# simulate coordinates
set.seed(22)
coords1 = cbind(runif(n[1],0,1), runif(n[1],0,1))
coords2 = cbind(runif(n[2],0,1), runif(n[2],0,1))
coords=rbind(coords1,coords2)

# Get U and V

# correlation function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}

# Definitions of matrix A, T and R
a=2
A <- matrix(0,q,q) 
A[1,1]=1
A[2,1]=0
A[2,2]=1

TT <- A%*%t(A)
T1=A[,1]%*%t(A[,1])
T2=A[,2]%*%t(A[,2])
#T1+T2==TT check
d=as.matrix(dist(coords))
R1=sigma1*apply(d,2,exp_corr,range=range1)
R2=sigma2*apply(d,2,exp_corr,range=range2)
RT1=kronecker(R1, T1)
RT2=kronecker(R2, T2)
R=RT1+RT2

# simulate Gaussian process
W=mvrnorm(n = 1, rep(0,80), R)
u1=W[c(1:n[1])]
u2=W[c((2*n[1]+1):((2*n[1])+n[2]))]
v2=W[c((2*n[1]+n[2]+1):length(W))]

# simulate response Y
Y1=u1+rnorm(n[1],tau1)
Y2=u2+v2+rnorm(n[2],tau2)

## to get the time series I should vary sigma1 sigma2 and a and loop. 



