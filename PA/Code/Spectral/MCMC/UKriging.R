# try prediction block
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


# just check prediction (not inference)
#CV 
K=5
npreds=50
U1.pred= array(0,dim=c(nrow(Y1),ncol(Y1),npreds))
RMSE=matrix(0,npreds,dim(Y1)[1])
COV=numeric(dim(Y1)[1])
U <- rbind(u1, u2)

for (i in 1:(dim(Y1)[1]/K))
{
  #set train and set sets of EPA data
  which.test1=seq(1:dim(Y1)[1])[(K*(i-1)+1):(K*i)]
  train.sites1=seq(1:dim(Y1)[1])[-which.test1]
  train.sites2 = c(train.sites1, seq((dim(Y1)[1]+1),(dim(Y1)[1])+dim(Y2)[1]))
  # Train-test split of locations
  Y11 = u1[train.sites1,]
  s11=s1[train.sites1,]
  sp1=s1[which.test1,]
  st1 = s1[train.sites1, ]
  st2 = s2
  # train set of PA (all)
  Y22 = u2
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
  all.d=as.matrix(dist(rbind(st1,st2,sp1,sp2)))
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

    # Ul is knwon
    if (F) {
      for (r in 1:nt)
      {
        Uls[,r]=t(chol(Mp00))%*%rnorm(n1+n2,0,sqrt(sigmau[r]))
      }
    }
    
    
    U1p = matrix(0,np1,nt)
    for (r in 1:nt)
    {
      Au=AA%*%as.matrix(U[train.sites2,r])
      sigmaB=sigmau[r]*B
      Ul.pred=rmvnorm(1,mean=Au,sigma=sigmaB)
      U1p[,r]=Ul.pred[(1:np1)]
    }
    
    U1.pred[which.test1,,j] <- beta.1+U1p
  }
  
  
  # Check prediction 
  
  for (k in 1:npreds)
  {
    a1=U1.pred[which.test1,,k]
    b1=u1[which.test1,]
    #RMSE
    RMSE[k,(K*(i-1)+1):(K*i)]=sqrt(apply((a1-b1)^2,1,sum,na.rm=TRUE)/dim(a1)[2])
  }
  
  m=1
  # Estimate coverage
  for(k in which.test1)
  {
    samps=U1.pred[m,,]
    b=u1[k,]
    
    a=apply(samps,1,quantile,c(0.05,0.95))
    cc=numeric(nt)
    for (j in 1:nt)
    {
      cc[j]=between(b[j],a[1,j],a[2,j])
    }
    
    COV[k]=mean(cc,na.rm=TRUE)
    m=m+1
  }
  
}

station=rep(1:(dim(Y1)[1]),each=npreds)
rmse=as.vector(RMSE[,1:(dim(Y1)[1])])
iterss=rep(1:npreds,(dim(Y1)[1]))

RMSES=data.frame(station,rmse,iterss)
RMSES$station=as.factor(RMSES$station)

station=seq(1:(dim(Y1)[1]))
cov=COV

COV.pr=data.frame(COV,station)
COV.pr$station=as.factor(COV.pr$station)


#RMSE
ggplot(RMSES)+geom_boxplot(aes(y=rmse,x=station))+theme_bw()
ggplot(RMSES %>% filter(station%in%seq(1,40)))+geom_boxplot(aes(y=rmse,x=station))+theme_bw()

p1=ggplot(RMSES)+geom_histogram(aes(x=rmse))+theme_bw()+ labs(x = "",y='')
p2=ggplot(RMSES)+geom_boxplot(aes(y=rmse))+theme_bw()
figureRMSE <- grid.arrange(p1, p2,ncol=2,top=textGrob("RMSE"))
