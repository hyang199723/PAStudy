
library(fields) 
library(geoR)
library(truncnorm)
source('ExtraFunctions.R')
source('LMC_function.R')

                ########################
                #### Simulated data ####
                ########################
#source('simAllTS.R')

                  ####################
                  #### Real Data #####
                  ####################
PA_data <- read.csv("Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv")
FRM_data <- read.csv("Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv")
# Convert timestamp
PA_data$Timestamp <- as.POSIXct(PA_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
FRM_data$Timestamp <- as.POSIXct(FRM_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
# No PA 
start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-01 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7

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

#exit2=Compact.LMC_fit(Y1,Y2, s1,s2)
#2042.504--> simudata 
exit1=LMC_fit(Y1,Y2, s1,s2)
#2389.701

# Plot results
par(mfrow=c(2,1))
plot(exit1$rangeU,type='l')
plot(exit2$rangeU,type='l')
abline(h=rangeu,col='red')  

par(mfrow=c(2,1))
plot(exit2$rangeV,type='l')
plot(exit1$rangeV,type='l')
abline(h=rangev,col='red')  

plot(exit1$tau1,type='l')
plot(exit2$tau1,type='l')
abline(h=tau1,col='red')  

plot(exit1$tau2,type='l')
plot(exit2$tau2,type='l')
abline(h=tau2,col='red')  

par(mfrow=c(3,4))
for (i in 1:dim(Y1)[2])
{
  plot(exit1$A[i,],type='l')
  #abline(h=al[i],col='red')  
}

par(mfrow=c(3,4))
for (i in 1:dim(Y1)[2])
{
  plot(exit1$sigmaU[i,],type='l')
  #abline(h=sigmau[i],col='red')  
}


par(mfrow=c(3,4))
for (i in 1:dim(Y1)[2])
{
  plot(exit2$sigmaV[i,],type='l')
  abline(h=sigmav[i],col='red')  
}

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



