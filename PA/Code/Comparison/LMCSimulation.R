rm(list  = ls())
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Comparison/")
source('LMC.R')
#### simulation of the time series process
library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)
library(mvtnorm)
library(truncnorm)
library(viridis)

# correlation function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}


# Constant correlation across all frequencies.
# Get strong correlation first.

###### set some parameters ########
set.seed(123)
# 350 + 100 Type1 and Type2 data
# 250 Type1 for training, 10 for testing
a1 = 100
a2 = 70
n=c(a1,a2) # number of locations
# Randomly select 100 testing locations
vld = rbinom(a1, 1, 0.302)

nt=15 # total time steps

tau1=6^2 # error variance1
tau2=8^2  # error variance2
set.seed(99)
#al=runif(nt,min = 5,max=6)
# Change al from uniform sequence to decreasing sequence
# from 10 - 1 / 10
al = seq(from = 0, to = 0, length = t) / 25 # /10: 20% - 80%;  /3: 40% - 98%; /25: 5% - 47%

# correlation parameters
set.seed(88)
sigmau=seq(from=10,to=1,length=nt)+rtruncnorm(nt,a=0,sd=.2)
set.seed(564)
sigmav=seq(from=5,to=.1,length=nt)+rtruncnorm(nt,a=0,sd=.2)
# same range for all freq
rangeu=3
lrangeu=log(rangeu)
rangev=2
lrangev=log(rangev)



###### simulate coordinates #######

set.seed(1)
leng = 15
coords1 = cbind(runif(n[1],0,leng), runif(n[1],0,leng))
set.seed(28)
coords2 = cbind(runif(n[2],0,leng), runif(n[2],0,leng))
coords=rbind(coords1,coords2)

## mean
lon1 = coords1[, 1]; lat1 = coords1[, 2]
lon1s = lon1 * lon1; lat1s = lat1 * lat1; lonlat1 = lon1 * lat1
X1 = cbind(rep(1, a1), lon1, lat1, lon1s, lat1s, lonlat1, rep(0, a1))

lon2 = coords2[, 1]; lat2 = coords2[, 2]
lon2s = lon2 * lon2; lat2s = lat2 * lat2; lonlat2 = lon2 * lat2
X2 = cbind(rep(1, a2), lon2, lat2, lon2s, lat2s, lonlat2, rep(1, a2))

set.seed(123)
beta = rnorm(7, 1, 2) / 100
beta[7] = sqrt(tau1) / 2

beta.1 = X1 %*% beta
beta.2 = X2 %*% beta

######## Get U and V ##########

du12=as.matrix(dist(coords)) #distance matrix U

u=matrix(NA,ncol=nt,nrow=sum(n))
u1=matrix(NA,ncol=nt,nrow=n[1])
u2=matrix(NA,ncol=nt,nrow=n[2])
v2=matrix(NA,ncol=nt,nrow=n[2])
dv2=as.matrix(dist(coords2)) # distance matrix v

M=exp_corr(du12, range=rangeu)
Sigmav22=exp_corr(dv2, range = rangev)
for (t in 1:nt)
{
  #u
  #M=exp_corr(du12,range=rangeu[t])
  u[,t]=t(chol(M))%*%rnorm(sum(n),0,sqrt(sigmau[t]))
  u1[,t]=u[(1:n[1]),t]
  u2[,t]=u[(n[1]+1):(sum(n)),t]
  
  #v
  #Sigmav22=exp_corr(dv2,range = rangev[t])
  v2[,t]=t(chol(Sigmav22))%*%rnorm(n[2],0,sqrt(sigmav[t]))
  
}

####### simulate response Y ############

#spectral
Z1sp=matrix(NA,ncol=nt,nrow=n[1])
Z2sp=matrix(NA,ncol=nt,nrow=n[2])

for (t in 1:nt)
{
  Z1sp[,t]=u1[,t] #+ rnorm(n[1],0,sqrt(nt/2*tau1))
  Z2sp[,t]=al[t]*u2[,t]+v2[,t] #+ rnorm(n[2],0,sqrt(nt/2*tau2))
}

#time domain
Y1=matrix(NA,ncol=nt,nrow=n[1])
Y2=matrix(NA,ncol=nt,nrow=n[2])
for(i in 1:n[1]){Y1[i,] <- fft_real(Z1sp[i,],inverse=TRUE)+beta.1[i]+rnorm(nt,0,sqrt(tau1))}
for(i in 1:n[2]){Y2[i,] <- fft_real(Z2sp[i,],inverse=TRUE)+beta.2[i]+rnorm(nt,0,sqrt(tau2))}

Y1_sum = cbind(Y1, vld)
Y1_train = subset(Y1, vld == 0)
Y1_test = subset(Y1, vld == 1)

c1 = cbind(coords1, vld)
coords1_train = subset(coords1, vld == 0)
coords1_test = subset(coords1, vld == 1)

cor = al*sigmau/sqrt((sigmau + tau1) * (al*al*sigmau + sigmav + tau2)) 
plot(x = 1:nt, y = cor, main = "correlation in spectral domain")
# Look at spatial correlation
plot(M[1, ], main = "spatial correlation, site 1")

# Save file
save(list = c("Y1_sum", "Y1_train", "Y1_test", "Y2", "coords1", "coords1_train", "coords1_test", "coords2"), file = "comparison.RData")

## Variance by time
v1 = rep(0, nt)
Y = rbind(Y1, Y2)
for (i in 1:nt) {
  v1[i] = var(Y[, i])
}
plot(v1, main = "Variance of LMC generated data")
