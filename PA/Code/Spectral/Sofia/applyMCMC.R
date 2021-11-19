# Script to apply LMC functions to real data
# Last update: 11/19/2021
rm(list=ls())
library(fields) 
library(geoR)
library(truncnorm)
library(tidyverse)
library(mvtnorm)
setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Spectral/ExtraFunctions.R')
source('Code/Spectral/LMC_function.R')

PA_data <- read.csv("Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv")
FRM_data <- read.csv("Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv")

# Convert timestamp
PA_data$Timestamp <- as.POSIXct(PA_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
FRM_data$Timestamp <- as.POSIXct(FRM_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
# No PA 
start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-02 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7

pa <- subset(PA_data, (Timestamp >= start) & (Timestamp <= end))
frm <- subset(FRM_data, (Timestamp >= start) & (Timestamp <= end))
# Get data to desired format
paTS <- pivot_wider(pa, names_from = Timestamp, values_from = PM25)
frmTS <- pivot_wider(frm, names_from = Timestamp, values_from = PM25)
# Record locations of PA and FRM stations
s1 <- as.matrix(frmTS[, 1:2])
s2 <- as.matrix(paTS[, 1:2])
# Get rid of the locations
paTS <- paTS[, -c(1:2)]
frmTS <- frmTS[, -c(1:2)]
Y1 = as.matrix(data.frame(frmTS))
colnames(Y1)=NULL
Y2 = as.matrix(data.frame(paTS))
colnames(Y2)=NULL
                        #####################
                        ### Fit the model ###
                        #####################
iters = 6000
thin = 1

time1 <- proc.time()[3]
exit1 = LMC_fit(Y1,Y2, s1,s2, iters=iters, thin=thin)
time2 <- proc.time()[3]
elap = time2 - time1

# Imputatoin quality
m1       <- is.na(Y1)
m2       <- is.na(Y2)

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

  
for (i in 1:thin) {
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
  val2=c(val2, as.vector(exit1$A[, , i]))
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

res1 = data.frame(val,param,nchain,iter,freq)
res2 = data.frame(val=val2,param=param2,nchain=nchain2,iter=iter2,freq=freq2)

res1$nchain=as.factor(res1$nchain)
res2$nchain=as.factor(res2$nchain)

#plot range
prang1=ggplot(res1 %>% filter(param=='rangeU')) + 
  geom_line(aes(x=iter,y=val,col=nchain)) +
  labs(title = 'range1') +
  theme_bw()
prang1
prang2=ggplot(res1 %>% filter(param=='rangeV')) +
  geom_line(aes(x=iter,y=val,col=nchain)) +
  labs(title = 'range2') +
  theme_bw()
prang2

ggsave('Figures/LMCConvergence/2Day/PostRange1.png', prang1)
ggsave('Figures/LMCConvergence/2Day/PostRange2.png', prang2)

ptau1=ggplot(res1 %>% filter(param=='tau1')) +
  geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5) +
  labs(title = 'tau1') +
  theme_bw()
ptau1
ptau2=ggplot(res1 %>% filter(param=='tau2')) +
  geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5) +
  labs(title = 'tau2') +
  theme_bw() +
  ylim(c(0,1))
ptau2
ggsave('Figures/LMCConvergence/2Day/PostTau1.png',ptau1)
ggsave('Figures/LMCConvergence/2Day/PostTau2.png',ptau2)

# plots by freq
freq1=1:round(dim(Y2)[2]/2)
freq2=(round(dim(Y2)[2]/2)+1):dim(Y2)[2]

pA_1=ggplot(res2 %>% filter(param=='A',freq %in% freq1)) +
  geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5) +
  facet_wrap(~freq,scales='free') +
  labs(title = 'A') +
  theme_bw()
pA_1
pA_2=ggplot(res2%>% filter(param=='A',freq %in% freq2)) +
  geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5) +
  facet_wrap(~freq,scales='free') +
  labs(title = 'A (cont.)') +
  theme_bw()
pA_2
ggsave('Figures/LMCConvergence/2Day/PostA_frq1:15.png',pA_1)
ggsave('Figures/LMCConvergence/2Day/PostA_frq16:30.png',pA_2)

psigmaU_1=ggplot(res2 %>% filter(param=='sigmaU',freq %in% freq1)) +
  geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5) +
  facet_wrap(~freq,scales='free') +
  labs(title = 'sigma1') +
  theme_bw()
psigmaU_1
psigmaU_2=ggplot(res2 %>% filter(param=='sigmaU',freq %in% freq2)) +
  geom_line(aes(x=iter,y=val,col=nchain),alpha=0.5) +
  facet_wrap(~freq,scales='free') +
  labs(title = 'sigma1 (cont.)') +
  theme_bw()
psigmaU_2

ggsave('Figures/LMCConvergence/2Day/PostSimaU_frq1:15.png',psigmaU_1)
ggsave('Figures/LMCConvergence/2Day/PostSimaU_frq16:30.png',psigmaU_2)

psigmaV_1=ggplot(res2 %>% filter(param=='sigmaV',freq %in% freq1)) +
  geom_line(aes(x=iter,y=val,col=nchain), alpha=0.5) +
  facet_wrap(~freq,scales='free') +
  labs(title = 'sigma2') +
  theme_bw()
psigmaV_1
psigmaV_2=ggplot(res2 %>% filter(param=='sigmaV',freq %in% freq2))+
  geom_line(aes(x=iter,y=val,col=nchain), alpha=0.5)+
  facet_wrap(~freq,scales='free')+
  labs(title = 'sigma2 (cont.)') +
  theme_bw()
psigmaV_2
ggsave('Figures/LMCConvergence/2Day/PostSigmaV_frq1:15.png',psigmaV_1)
ggsave('Figures/LMCConvergence/2Day/PostSigmaV_frq16:30.png',psigmaV_2)


########################################################
############# Cross validation analysis ################
########################################################

nsites.1=nrow(Y1)
nsites.2=nrow(Y2)
burn=3000
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


for (i in 1:K) {
  r.n = round(nsites.1 / K)
  which.test1=((i-1) * r.n + 1):(i * r.n)
  train.sites1=seq(1:nsites.1)[-which.test1]
  
  r.n = round(nsites.2 / K)
  which.test2 = ((i-1) * r.n + 1):(i * r.n)
  train.sites2 = seq(1:nsites.2)[-which.test2]

  Y11 = Y1
  Y11[which.test1,]=NA
  
  Y22 = Y2
  Y22[which.test2,]=NA
  
  exit=LMC_fit(Y11,Y22, s1,s2,iters=iters,burn=burn)
  
  for (j in burn:iters) {
    new.vals1=(exit$Y1.m[which.test1, , j, 1] - Y1[which.test1,])^2
    vals1=c(vals1,as.vector(new.vals1))
    times1=c(times1,rep(1:dim(new.vals1)[2], each=dim(new.vals1)[1]))
    sites1=c(sites1,rep(which.test1,dim(new.vals1)[2]))
  
    new.vals2=(exit$Y2.m[which.test2, , j, 1]-Y2[which.test2,])^2
    vals2=c(vals2,as.vector(new.vals2))
    times2=c(times2,rep(1:dim(new.vals2)[2],each=dim(new.vals2)[1]))
    sites2=c(sites2,rep(which.test2,dim(new.vals2)[2]))
  }
  
  posIter1=c(posIter1, rep(seq(burn:iters), prod(dim(new.vals1))))
  ks1=c(ks1,rep(i, each = prod(dim(new.vals1)) * (iters + 1 - burn)))
  
  posIter2=c(posIter2, rep(seq(burn:iters), prod(dim(new.vals2))))
  ks2=c(ks2,rep(i, each=prod(dim(new.vals2))*(iters+1-burn)))
}



cv1 = data.frame(ks1, posIter1, times1, sites1, vals1)
cv2 = data.frame(ks2, posIter2, times2, sites2, vals2)


RMSE1=cv1 %>% group_by(ks1,sites1,posIter1) %>% summarise(RMSE=sqrt((mean(vals1,na.rm=TRUE))))
RMSE1
RMSE2=cv2 %>% group_by(ks2,sites2,posIter2) %>% summarise(RMSE=sqrt((mean(vals2,na.rm=TRUE))))
RMSE2

RMSE1$sites1=as.factor(RMSE1$sites1)
RMSE2$sites2=as.factor(RMSE2$sites2)

plot1=ggplot(RMSE1)+geom_boxplot(aes(y=RMSE,x=sites1))+theme_bw()
plot1
plot2=ggplot(RMSE2)+geom_boxplot(aes(y=RMSE,x=sites2))+theme_bw()
plot2


ggsave("Figures/ImputationCV/2Day/RMSEVC1.pdf",plot=plot1)
ggsave("Figures/ImputationCV/2Day/RMSEVC2.pdf",plot=plot2)