# This script is for computationally intense tasks that runs on high performance machine
# All data and source files need to be stored under the same directory
# This includes the main LMC functions
# This script has two functions, LMC_fit and compact.LMC_fit
# Last update: 12/01/2021

# Script to apply LMC functions to real data
# Last update: 12/01/2021
rm(list=ls())
# Set working directory to current
setwd('/Users/hongjianyang/Research/PAStudy/PA/')

library(fields) 
library(glue)
library(viridis)
#library(geoR)
library(truncnorm)
library(tidyr)
library(mvtnorm)
library(ggplot2)
source('Code/Spectral/ExtraFunctions.R')
source('Code/Spectral/LMC_function.R')
frmTS <- read.csv('Data/Formatted_PA_FRM/missing_FRM.csv')
paTS <- read.csv('Data/Formatted_PA_FRM/missing_PA.csv')
frm.impute <- read.csv('Data/Formatted_PA_FRM/EPA_Imputed_2020.csv')
pa.impute <- read.csv('Data/Formatted_PA_FRM/PA_Imputed_2020.csv')

s1 <- as.matrix(frmTS[, 1:2])
s2 <- as.matrix(paTS[, 1:2])
# Get rid of the locations
frmTS <- frmTS[, -c(1:2)]
paTS <- paTS[, -c(1:2)]
Y1 = as.matrix(data.frame(frmTS))
colnames(Y1)=NULL
Y2 = as.matrix(data.frame(paTS))
colnames(Y2)=NULL

# Make one complete station 0
rowNum = which(rowSums(is.na(Y2)) == 0)[1] # Row Number is 7
raw <- Y2[rowNum, ]
Y2[rowNum, ] = NA


#####################
### Fit the model ###
#####################

#exit2=Compact.LMC_fit(Y1,Y2, s1,s2,iters=6000)
#2042.504--> simudata 
start = proc.time()[3]
iters = 6000
thin = 1
exit1 = LMC_fit(Y1, Y2, s1, s2, iters = iters, thin = thin)
end = proc.time()[3]
#2389.701

# Get Y1 and Y2
burnin = 2000
y1.raw <- exit1$Y1.m[,,burnin:iters,1]
y2.raw <- exit1$Y2.m[,,burnin:iters,1]

y1.complete = rowMeans(y1.raw, dims = 2)
y2.complete = rowMeans(y2.raw, dims = 2)

imp <- y2.complete[rowNum, ]
cor(imp, raw)

                    #################################
                    ### Analyse Exit of the model ###
                    #################################

# Create data frame
nchain=c()
param=c()
val=c()
iter=c()
freq=c()

nchain2=c()
param2=c()
val2=c()
iter2=c()
freq2=c()

for (i in 1:thin)
{
  #rangeU
  val=c(val,exit1$rangeU[,i])
  param=c(param,rep("rangeU",iters))
  nchain=c(nchain,rep(i,iters))
  iter=c(iter,1:iters)
  freq=c(freq,rep(NA,iters))
  #rangeV
  val=c(val,exit1$rangeV[,i])
  param=c(param,rep("rangeV",iters))
  nchain=c(nchain,rep(i,iters))
  iter=c(iter,1:iters)
  freq=c(freq,rep(NA,iters))
  #tau1 
  val=c(val,exit1$tau1[,i])
  param=c(param,rep("tau1",iters))
  nchain=c(nchain,rep(i,iters))
  iter=c(iter,1:iters)
  freq=c(freq,rep(NA,iters))
  #tau2
  val=c(val,exit1$tau2[,i])
  param=c(param,rep("tau2",iters))
  nchain=c(nchain,rep(i,iters))
  iter=c(iter,1:iters)
  freq=c(freq,rep(NA,iters))
  #A
  val2=c(val2,as.vector(exit1$A[,,i]))
  param2=c(param2,rep('A',ncol(Y2)*iters))
  nchain2=c(nchain2,rep(i,ncol(Y2)*iters))
  iter2=c(iter2,rep(1:iters,each=ncol(Y2)))
  freq2=c(freq2,rep(1:ncol(Y2),iters))
  #sigmaU
  val2=c(val2,as.vector(exit1$sigmaU[,,i]))
  param2=c(param2,rep('sigmaU',ncol(Y2)*iters))
  nchain2=c(nchain2,rep(i,ncol(Y2)*iters))
  iter2=c(iter2,rep(1:iters,each=ncol(Y2)))
  freq2=c(freq2,rep(1:ncol(Y2),iters))
  #sigmaV
  val2=c(val2,as.vector(exit1$sigmaV[,,i]))
  param2=c(param2,rep('sigmaV',ncol(Y2)*iters))
  nchain2=c(nchain2,rep(i,ncol(Y2)*iters))
  iter2=c(iter2,rep(1:iters,each=ncol(Y2)))
  freq2=c(freq2,rep(1:ncol(Y2),iters))
  
}

res1=data.frame(val,param,nchain,iter,freq)
res2=data.frame(val=val2,param=param2,nchain=nchain2,iter=iter2,freq=freq2)

res1$nchain=as.factor(res1$nchain)
res2$nchain=as.factor(res2$nchain)

#plot range
prang1=ggplot(res1 %>% filter(param=='rangeU'))+geom_line(aes(x=iter,y=val,col=nchain))+theme_bw()
prang2=ggplot(res1 %>% filter(param=='rangeV'))+geom_line(aes(x=iter,y=val,col=nchain))+theme_bw()

# ggsave('PostRange1.png',prang1)
# ggsave('PostRange2.png',prang2)

ptau1=ggplot(res1 %>% filter(param=='tau1'))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+theme_bw()
ptau2=ggplot(res1 %>% filter(param=='tau2'))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+theme_bw()+ylim(c(0,1))

# ggsave('PostTau1.png',ptau1)
# ggsave('PostTau2.png',ptau2)

# plots by freq
freq1=1:round(dim(Y2)[2]/2)
freq2=(round(dim(Y2)[2]/2)+1):dim(Y2)[2]

pA_1=ggplot(res2%>% filter(param=='A',freq %in%freq1))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+facet_wrap(~freq,scales='free')+theme_bw()
pA_2=ggplot(res2%>% filter(param=='A',freq %in%freq2))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+facet_wrap(~freq,scales='free')+theme_bw()

# ggsave('PostA_frq1:15.png',pA_1)
# ggsave('PostA_frq16:30.png',pA_2)

psigmaU_1=ggplot(res2%>% filter(param=='sigmaU',freq %in%freq1))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+facet_wrap(~freq,scales='free')+theme_bw()
psigmaU_2=ggplot(res2%>% filter(param=='sigmaU',freq %in%freq2))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+facet_wrap(~freq,scales='free')+theme_bw()

# ggsave('PostSimaU_frq1:15.png',psigmaU_1)
# ggsave('PostSimaU_frq16:30.png',psigmaU_2)

psigmaV_1=ggplot(res2%>% filter(param=='sigmaV',freq %in%freq1))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+facet_wrap(~freq,scales='free')+theme_bw()
psigmaV_2=ggplot(res2%>% filter(param=='sigmaV',freq %in%freq2))+geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5)+facet_wrap(~freq,scales='free')+theme_bw()

# ggsave('PostSigmaV_frq1:15.png',psigmaV_1)
# ggsave('PostSigmaV_frq16:30.png',psigmaV_2)


########################################################
############# Cross validation analysis ################
########################################################

nsites.1=nrow(Y1)
nsites.2=nrow(Y2)
iters=3000
burn=1000
K=5

ks1=c()
sites1=c()
times1=c()
vals1=c()
posIter1=c()

ks2=c()
sites2=c()
times2=c()
vals2=c()
posIter2=c()


for (i in 1:K)
{
  r.n=round(nsites.1/K)
  which.test1=((i-1)*r.n+1):(i*r.n)
  train.sites1=seq(1:nsites.1)[-which.test1]
  
  r.n=round(nsites.2/K)
  which.test2=((i-1)*r.n+1):(i*r.n)
  train.sites2=seq(1:nsites.2)[-which.test2]

  Y11 = Y1
  Y11[which.test1,]=NA
  
  Y22 = Y2
  Y22[which.test2,]=NA
  
  exit=LMC_fit(Y11,Y22, s1,s2,iters=iters,burn=burn)
  
  for (j in burn:iters)
  {
    
    new.vals1=(exit$Y1.m[which.test1,,j]-Y1[which.test1,])^2
    vals1=c(vals1,as.vector(new.vals1))
    times1=c(times1,rep(1:dim(new.vals)[2],each=dim(new.vals)[1]))
    sites1=c(sites1,rep(which.test1,dim(new.vals)[2]))
  
    new.vals2=(exit$Y2.m[which.test2,,j]-Y2[which.test2,])^2
    vals2=c(vals2,as.vector(new.vals2))
    times2=c(times2,rep(1:dim(new.vals2)[2],each=dim(new.vals2)[1]))
    sites2=c(sites2,rep(which.test2,dim(new.vals2)[2]))
    
  }
  posIter1=c(posIter1,rep(seq(burn:iters),prod(dim(new.vals1))))
  ks1=c(ks1,rep(i,each=prod(dim(new.vals1))*(iters+1-burn)))
  
  posIter2=c(posIter2,rep(seq(burn:iters),prod(dim(new.vals2))))
  ks2=c(ks2,rep(i,each=prod(dim(new.vals2))*(iters+1-burn)))
}



cv1=data.frame(ks1,posIter1,times1,sites1,vals1)
cv2=data.frame(ks2,posIter2,times2,sites2,vals2)

library(tidyverse)

RMSE1=cv1 %>% group_by(ks1,sites1,posIter1) %>% summarise(RMSE=sqrt((mean(vals1,na.rm=TRUE))))
RMSE2=cv2 %>% group_by(ks2,sites2,posIter2) %>% summarise(RMSE=sqrt((mean(vals2,na.rm=TRUE))))

RMSE1$sites1=as.factor(RMSE1$sites1)
RMSE2$sites2=as.factor(RMSE2$sites2)

plot1=ggplot(RMSE1)+geom_boxplot(aes(y=RMSE,x=sites1))+theme_bw()
plot2=ggplot(RMSE2)+geom_boxplot(aes(y=RMSE,x=sites2))+theme_bw()


# ggsave("RMSEVC1.pdf",plot=plot1)
# ggsave("RMSEVC2.pdf",plot=plot2)



