# try prediction block
rm(list = ls())
library(fields) 
library(geoR)
library(truncnorm)
library(tidyverse)
library(mvtnorm)
source('ExtraFunctions.R')
source('LMC_function.R')

#simulate data
source('simAllTS.R') # load your data here 


RMSE=matrix(NA,ncol=dim(Y1)[1],nrow=2000)
COV=numeric(dim(Y1)[1])
iters=3000
burn=1000


#CV 
K=10 

for (i in 1:(dim(Y1)[1]/K))
{
  #set train and set sets of EPA data
  which.test1=seq(1:dim(Y1)[1])[(K*(i-1)+1):(K*i)]
  train.sites1=seq(1:dim(Y1)[1])[-which.test1]
  
  Y11 = Y1[train.sites1,]
  s11=s1[train.sites1,]
  sp1=s1[which.test1,]
  
  # train set of PA (all)
  Y22 = Y2
  s22=s2
  sp2=NULL # no test set
  
  # fit model and predict
  exit=LMC_fit(Y1=Y11,Y2=Y22, s1=s11,s2=s22,sp1=sp1,sp2=sp2,
               mean_range=0, sd_range=1, mean_var=0, sd_var=1, mean_rho=0,
               sd_rho=10, iters=iters, burn=burn, thin=1, update=10)
  
  # Check prediction 
  
  for (k in 2:(iters-burn))
  {
    a1=exit$Y1p[,,k]  
    b1=Y1[which.test1,]
    #RMSE
    RMSE[k,(K*(i-1)+1):(K*i)]=sqrt(apply((a1-b1)^2,1,sum,na.rm=TRUE)/dim(a1)[2])
  }
  
  
  # Estimate coverage
  for(k in (K*(i-1)+1):(K*i))
  {
    n.site=k
    samps=exit$Y1p[n.site,,]
    b=Y1[n.site,]
    
    a=apply(samps,1,quantile,c(0.05,0.95))
    cc=numeric(nt)
    for (j in 1:nt)
    {
      cc[j]=between(b[j],a[1,j],a[2,j])
    }
    
    COV[n.site]=mean(cc,na.rm=TRUE)
  }
  #print(i)  
}

# get data frame to plot
station=rep(1:40,each=2000)
rmse=as.vector(RMSE[,1:40])
iterss=rep(1:2000,40)

RMSES=data.frame(station,rmse,iterss)
RMSES$station=as.factor(RMSES$station)

ggplot(RMSES)+geom_boxplot(aes(y=rmse,x=station))+theme_bw()



