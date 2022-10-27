# rm(list  = ls())
# setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Simulation/")
source('LMC_Sim.R')
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
a1 = 80
a2 = 300
n=c(a1,a2) # number of locations
# Randomly select 100 testing locations
vld = rbinom(a1, 1, 0.302)

nt=91 # total time steps
ntot=nt*(a1+a2)

tau1=6^2 # error variance1
tau2=8^2  # error variance2
set.seed(99)
#al=runif(nt,min = 5,max=6)
# Change al from uniform sequence to decreasing sequence
# from 10 - 1 / 10
al = seq(from = 0, to = 0, length = nt) / 25 # /10: 20% - 80%;  /3: 40% - 98%; /25: 5% - 47%

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

## Mean: create covariates
X=array(NA,dim=c(a1+a2,nt,6))
X[,,1]=matrix(rep(1,ntot),ncol = nt)
X[,,2]=matrix(rnorm(ntot,mean=47,sd=sqrt(316)),ncol = nt)
X[,,3]=matrix(rnorm(ntot,mean=66,sd=sqrt(684)),ncol = nt)
X[,,4]=matrix(rbinom(ntot,size=1,prob=.4),ncol = nt)
X[,,5]=matrix(rbinom(ntot,size=1,prob=.3),ncol = nt)
X[,,6]=matrix(rbinom(ntot,size=1,prob=.1),ncol = nt)

X1=X[1:a1,,]
X2=X[(a1+1):(a1+a2),,]

set.seed(123)
betau = rnorm(6, 1, 2) / 100
betav = rnorm(6, .5, 2) / 100
#beta[7] = sqrt(tau1) / 2

# #X\beta
Xbu=array(NA,dim=c(a1+a2,nt,6))
Xbv=array(NA,dim=c(a2,nt,6))

for (i in 1:6)
{
  Xbu[,,i]=X[,,i]*betau[i]
  Xbv[,,i]=X2[,,i]*betav[i]
}

# Get DFT
Xbeta.uT=array(NA,dim=c(a1+a2,nt,6))
Xbeta.vT=array(NA,dim=c(a2,nt,6))

for (i in 1:6)
{Xbeta.uT[,,i]=t(apply(Xbu[,,i],1,fft_real))
Xbeta.vT[,,i]=t(apply(Xbv[,,i],1,fft_real))}

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
  u1[,t]=u[(1:n[1]),t]+apply(Xbeta.uT[(1:a1),t,],1,sum)
  u2[,t]=u[(n[1]+1):(sum(n)),t]+apply(Xbeta.uT[(a1+1):(a1+a2),t,],1,sum)
  
  #v
  #Sigmav22=exp_corr(dv2,range = rangev[t])
  v2[,t]=t(chol(Sigmav22))%*%rnorm(n[2],0,sqrt(sigmav[t]))+apply(Xbeta.vT[,t,],1,sum)
  
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
for(i in 1:n[1]){Y1[i,] <- fft_real(Z1sp[i,],inverse=TRUE)+rnorm(nt,0,sqrt(tau1))}
for(i in 1:n[2]){Y2[i,] <- fft_real(Z2sp[i,],inverse=TRUE)+rnorm(nt,0,sqrt(tau2))}

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


## Some plots

valuesY=c(as.vector(Y1),Y1=as.vector(Y2))
type=c(rep('Type1',length(Y1)),rep('Type2',length(Y2)))
xcord=c(rep(coords1[,1],nt),rep(coords2[,1],nt))
ycord=c(rep(coords1[,2],nt),rep(coords2[,2],nt))
times=c(rep(seq(1,nt),each=a1),rep(seq(1,nt),each=a2))

df=data.frame(valuesY,type,xcord,ycord,times)


ggplot(df %>% filter(times==1))+geom_point(aes(x=xcord,y=ycord,col=valuesY))+
  theme_bw()+facet_grid(~type)

ggplot(df %>% filter(times==4))+geom_point(aes(x=xcord,y=ycord,col=valuesY))+
  theme_bw()+facet_grid(~type)


