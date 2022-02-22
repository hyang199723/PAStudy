# This script is for computationally intense tasks that runs on high performance machine
# All data and source files need to be stored under the same directory
# This includes the main LMC functions
# This script has two functions, LMC_fit and compact.LMC_fit
# Last update: 12/01/2021

# Script to apply LMC functions to real data
# Last update: 12/01/2021
rm(list=ls())
# Set working directory to current
library(rstudioapi)
curpath = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curpath)

library(fields) 
library(glue)
library(viridis)
#library(geoR)
library(truncnorm)
library(tidyr)
library(mvtnorm)
library(ggplot2)
source('ExtraFunctions.R')
source('LMC_function.R')
OR = as.POSIXct('1970-01-01', tz = 'UTC')

########################
#### Simulated data ####
########################
#source('simAllTS.R')

####################
#### Real Data #####
####################
PA_raw <- read.csv("PA_2020_Hourly_Formatted.csv")
FRM_raw <- read.csv("FRM_2020_Hourly_Formatted.csv")
PA_raw <- PA_raw[order(PA_raw$Timestamp), ]
FRM_raw <- FRM_raw[order(FRM_raw$Timestamp), ]
PA_data <- PA_raw
FRM_data <- FRM_raw
# Convert timestamp
# Note: We need to convert time to EST to avoid NA time issue
# For PA data, there is no missingness after conversion
# However, for FRM data, there are 20 missing values
PA_data$Timestamp <- as.POSIXct(PA_raw$Timestamp, tz = 'EST', format = "%Y-%m-%d %H:%M:%OS")
FRM_data$Timestamp <- as.POSIXct(FRM_raw$Timestamp, tz= 'EST', format = "%Y-%m-%d %H:%M:%OS")

# Get rid of missingness in FRM
FRM_data <- FRM_data[!is.na(FRM_data$Timestamp), ]
# Complete Purple Air timestamps
time <- unique(PA_data$Timestamp)
int.seq <- as.numeric(time) / 3600
int.start = int.seq[1]
int.end = int.seq[length(int.seq)]
complete = int.start:int.end
pa.missing = complete[!(complete %in% int.seq)]
time.missing = as.POSIXct(pa.missing * 3600, origin = OR)
lon = PA_data[1,1]
lat = PA_data[1,2]
df <- data.frame(Lon = lon, Lat = lat, Timestamp = time.missing, PM25 = NA)
PA.complete = rbind(PA_data, df)

# Complete FRM data
time <- unique(FRM_data$Timestamp)
int.seq <- as.numeric(time) / 3600
int.start = int.seq[1]
int.end = int.seq[length(int.seq)]
complete = int.start:int.end
frm.missing = complete[!(complete %in% int.seq)]
time.missing = as.POSIXct(frm.missing * 3600, origin = OR)
lon = FRM_data[1,1]
lat = FRM_data[1,2]
df <- data.frame(Lon = lon, Lat = lat, Timestamp = rep(time.missing, each = length(lon)), PM25 = NA)
FRM.complete = rbind(FRM_data, df)

# Check NAs in PA and FRM data
sum(is.na(PA.complete$Timestamp))
sum(is.na(FRM.complete$Timestamp))

start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-03 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7
interval = (as.numeric(end) - as.numeric(start)) / 3600
print(interval + 1)

pa <- subset(PA.complete, (Timestamp >= start) & (Timestamp <= end))
frm <- subset(FRM.complete, (Timestamp >= start) & (Timestamp <= end))
# Get data to desired format
frmTS <- pivot_wider(frm, names_from = Timestamp, values_from = PM25)
paTS <- pivot_wider(pa, names_from = Timestamp, values_from = PM25)

dim(frmTS)
dim(paTS)
# Record locations of PA and FRM stations
s1 <- as.matrix(frmTS[, 1:2])
s2 <- as.matrix(paTS[, 1:2])
# Get rid of the locations
frmTS <- frmTS[, -c(1:2)]
paTS <- paTS[, -c(1:2)]
Y1 = as.matrix(data.frame(frmTS))
colnames(Y1)=NULL
Y2 = as.matrix(data.frame(paTS))
colnames(Y2)=NULL
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
#write.csv(y1.complete, 'EPA_Imputed_2020.csv')
#write.csv(y2.complete, 'PA_Imputed_2020.csv')

y1.miss <- is.na(Y1)
y2.miss <- is.na(Y2)
## Make plots of missing values for a few stations
# First, make plots for EPA stations

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
thin = 1
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






if (F) {
  mean_range=0
  sd_range=1
  mean_var=0
  sd_var=1
  mean_rho=0
  sd_rho=10
  iters=6000
  burn=1000
  thin=1
  update=10
  
  
  Y1 = as.matrix(Y1)
  Y2 = as.matrix(Y2)
  
  n1       <- nrow(Y1)
  n2       <- nrow(Y2)
  nt       <- ncol(Y1)
  m1       <- is.na(Y1)
  m2       <- is.na(Y2)
  d        <- as.matrix(dist(rbind(s1,s2)))
  dv2 = as.matrix(dist(s2))
  const    <- nt/2
  
  # Keep track of stuff
  
  # All parameters: "rangeu","rangev","sigmaU","sigmaV","taue1","taue2","A",'beta1','beta2'
  
  # Parameters that are different for each spectrum: sigmaU, sigmaV, A, beta1, beta2
  # Dimension of data structure: nt * iterations
  keep.sigmaU = array(0, dim = c(nt, iters,thin))
  keep.sigmaV = array(0, dim = c(nt, iters,thin))
  keep.A = array(0, dim = c(nt, iters,thin))
  
  # Parameters that stay the same: rangeU, rangeV, taus1, taus2
  # Dimension of data structure: 1 * iterations
  keep.rangeU = matrix(0, iters,thin)
  keep.rangeV = matrix(0, iters,thin)
  keep.taus1 = matrix(0, iters,thin)
  keep.taus2 = matrix(0, iters,thin)
  keep.taue1 = matrix(0, iters,thin)
  keep.taue2 = matrix(0, iters,thin)
  # Vector parameters; Dimension of data structure: #stations * nt * iters
  # keep.u1= array(0,dim=c(n1,nt,iters))
  # keep.u2= array(0,dim=c(n2,nt,iters))
  # keep.v2= array(0,dim=c(n2,nt,iters))
  keep.Y1.M= array(0,dim=c(n1,nt,iters,thin))
  keep.Y2.M= array(0,dim=c(n2,nt,iters,thin))
  
  
  ## get info for the ranges priors 
  
  priorR_mn1 <- log(max(d)) - 1.5
  priorR_sd1 <- 1
  
  priorR_mn2 <- log(max(dv2)) - 1.5
  priorR_sd2 <- 1
  
  
  
  # start MCMC
  start = proc.time()[3]
  for(ttt in 1:thin){
    # Mean imputation
    beta1 <- mean(Y1,na.rm=TRUE)
    beta2 <- mean(Y2,na.rm=TRUE)
    Y1[m1] <- rnorm(sum(m1), beta1, 3) # Random variance to make it work
    Y2[m2] <- rnorm(sum(m2), beta2, 3) 
    # create initial vectors 
    U1     <- matrix(0,n1,nt)
    U2     <- matrix(0,n2,nt)
    V2     <- matrix(0,n2,nt)
    Z1     <- matrix(0,n1,nt)
    Z2     <- matrix(0,n2,nt)
    Ys1    <- matrix(0,n1,nt)
    Ys2    <- matrix(0,n2,nt)
    # Initial parameter values
    rangeU  <- 0.2
    rangeV   <- 0.2
    lrangeU  <- log(rangeU)
    lrangeV   <-log(rangeV)
    taue1  <- var(as.vector(Y1),na.rm=TRUE)
    taue2  <- var(as.vector(Y2),na.rm=TRUE)
    sigmaU   <- rep(taue1,nt)
    sigmaV   <- rep(taue1,nt)
    A      <- rep(0,nt)
    
    
    for(iter in 1:iters){
      ##############################################:
      ####     Transform to spatial land       #####:
      ##############################################:
      
      for(i in 1:n1){Z1[i,] <- fft_real(U1[i,],inverse=TRUE)}
      for(i in 1:n2){Z2[i,] <- fft_real(A*U2[i,]+V2[i,],inverse=TRUE)}
      
      ##############################################:
      ####  IMPUTE MISSING DATA (real space)   #####:
      ##############################################:
      
      Y1[m1] <- rnorm(sum(m1),beta1+Z1[m1],sqrt(taue1)) 
      Y2[m2] <- rnorm(sum(m2),beta2+Z2[m2],sqrt(taue2)) 
      
      ##############################################:
      ####      MEAN PARAMETERS (real space)   #####:
      ##############################################:
      
      # full conditional for beta1 and beta2 ...
      #VVV   <- (taue1*n1*nt + 0.01)
      #MMM   <- taue1*sum(Y1-Z1) 
      #beta1 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      #beta1=beta.1
      
      # VVV   <- (taue2*n2*nt + 0.01)
      # MMM   <- taue2*sum(Y2-Z2) 
      #beta2 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      #beta2=beta.2
      
      ##############################################:
      ####     ERROR VARIANCES (real space)    #####:
      ##############################################:
      # full conditionals for taue1 and taue2 
      taue1 <- 1/rgamma(1,n1*nt/2+1,sum((Y1-beta1-Z1)^2)/2+1)
      taue2 <- 1/rgamma(1,n2*nt/2+1,sum((Y2-beta2-Z2)^2)/2+1)
      
      ##############################################:
      ####     Transform to spectral land      #####:
      ##############################################:
      
      for(i in 1:n1){Ys1[i,] <- fft_real(as.numeric(Y1[i,] - beta1))}
      for(i in 1:n2){Ys2[i,] <- fft_real(as.numeric(Y2[i,] - beta2))}
      taus1 <- nt / 2 *taue1
      taus2 <- nt / 2 *taue2
      
      ##############################################:
      ####      LMC TERMS (spectral space)     #####:
      ##############################################:
      
      eigU = invU(d,n1,n2,rangeU)
      S11 = eigU$S11; S12 <- eigU$S12; S21 <- eigU$S21; S22 <- eigU$S22
      S1 = S11 - eigU$A12 %*% S21
      ES1=eigen(S1)
      S1_G   = ES1$vectors
      S1_D   = ES1$values
      S1inv = S1_G%*%diag(1/S1_D)%*%t(S1_G)
      A1 = S12 %*% solve(S22)
      S1invA1 = S1inv %*% A1
      # S1inv = G * diag(1/D) * G'
      # inverse = G * [1/taus1 + 1/sigmaU * diag(1/D)]^{-1} * G'
      
      S11inv=solve(S11)
      S2 = S22 - S21 %*% S11inv %*% S12
      A2 = S21 %*% S11inv
      ES2 = eigen(S2)
      S2_G = ES2$vectors
      S2_D = ES2$values
      S2inv = S2_G%*%diag(1/S2_D)%*%t(S2_G)
      S2invA2 = S2inv %*% A2
      
      # G: eigenvectors
      # D: eigenvalues
      # Q: inverse matrix
      eigV  = invV(dv2,rangeV)
      V_G = eigV$G
      V_D = eigV$D
      
      S2Vinv= invV(dv2,rangeV)$Q
      
      # Gibbs sampling
      for (r in 1:nt) # for each spectral 
      {
        # # Sample U1
        newDiag = 1/taus1 + 1/sigmaU[r] * 1/S1_D
        sigmaU1 <- S1_G %*% diag(1/newDiag) %*% t(S1_G)
        # S1_G
        meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1invA1 %*% U2[,r])
        U1[,r] = meanU1+S1_G%*%rnorm(n1, 0, 1/sqrt(newDiag))
        
        # # Sample U2
        newDiag = A[r]^2/taus2 + 1/sigmaU[r] * 1/S2_D
        sigmaU2 <- S2_G %*% diag(1/newDiag) %*% t(S2_G)
        
        meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 *A[r] * V2[,r] +
                                 1/sigmaU[r] * S2invA2 %*% U1[,r])
        U2[,r] = meanU2+S2_G%*%rnorm(n2,0,1/sqrt(newDiag))
        
        # # Sample V2
        newDiag = 1/sigmaV[r] * 1/V_D + 1/taus2
        sigmaV2 = V_G %*% diag(1/newDiag) %*% t(V_G)
        
        meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
        V2[,r] =meanV2+V_G%*%rnorm(n2,0,1/sqrt(newDiag))      
        
        # # Sample Al
        sigmaAl <- solve(1/taus2 * t(U2[,r]) %*% U2[,r] + 1/5)
        meanAl <- sigmaAl %*% (1/taus2 * t(Ys2[,r]) %*% U2[,r] - 1/taus2 * t(V2[,r]) %*% U2[,r] + 0.8/5)
        A[r] <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl, sd=sqrt(sigmaAl))
        
        # # Sample sig1
        U_sim <- as.vector(append(U1[,r], U2[,r]))
        a <- (n2+n1)/2 + 1
        b <- (t(U_sim) %*% eigU$Q %*% U_sim) / 2 + 1
        sigmaU[r] = 1/rgamma(1, a, b)
        
        # # # Sample sig2
        a <- n2/2 + 1
        b <- (t(V2[,r]) %*% eigV$Q %*% V2[,r]) / 2 + 1 
        sigmaV[r] = 1/rgamma(1, a, b)
        
      }
      
      ###############################################
      ##       Metropolis H: Range parameters      ##
      ###############################################
      # Sweep operation to get rid of variance
      Ru=sweep(rbind(U1,U2),2,FUN='/',sqrt(sigmaU))
      Rv=sweep(V2,2,FUN='/',sqrt(sigmaV))
      
      # range1
      Ms=exp_corr(d,range=exp(lrangeU))
      #curll = sum(apply(Ru,2,dmvnorm,mean=rep(0,n1+n2),sigma=Ms,log=TRUE))
      curll = sum(dmvnorm(t(Ru), mean=rep(0,n1+n2),sigma=Ms,log=TRUE))
      prior_curll=dnorm(lrangeU,mean=priorR_mn1,sd=priorR_sd1,log=TRUE)
      
      canrange1 = rnorm(1,lrangeU,0.1)
      canM = exp_corr(d,range=exp(canrange1))
      #canll = sum(apply(Ru,2,dmvnorm,mean=rep(0,n1+n2),sigma=canM,log=TRUE))
      canll = sum(dmvnorm(t(Ru), mean=rep(0,n1+n2),sigma=canM,log=TRUE))
      prior_canll=dnorm(canrange1,mean=priorR_mn1,sd=priorR_sd1,log=TRUE)
      
      print(iter)
      MH1 <- canll-curll+prior_canll-prior_curll+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
      
      
      
      if (log(runif(1))<MH1)
      {
        lrangeU=canrange1
      }
      
      # range2
      Ss=exp_corr(dv2,range = exp(lrangeV))
      #curll2 = sum(apply(Rv,2,dmvnorm,mean=rep(0,n2),sigma=Ss,log=TRUE))
      curll2 = sum(dmvnorm(t(Rv), mean = rep(0,n2),sigma=Ss,log=TRUE ))
      prior_curll2=dnorm(lrangeV,mean=priorR_mn2,sd=priorR_sd2,log=TRUE)
      
      canrange2 = rnorm(1,lrangeV,0.1)
      canS = exp_corr(dv2,range=exp(canrange2))
      #canll2 = dmvnorm(Rv,rep(0,n2),canS,log=TRUE)
      canll2 = sum(dmvnorm(t(Rv), mean = rep(0,n2),sigma=canS,log=TRUE ))
      prior_canll2=dnorm(canrange2,mean=priorR_mn2,sd=priorR_sd2,log=TRUE)
      
      MH2 <- canll2-curll2+prior_canll2-prior_curll2+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
      
      if (log(runif(1))<MH2)
      {
        lrangeV=canrange2
      }
      rangeU = exp(lrangeU)
      rangeV = exp(lrangeV)
      
      
      
      ##############################################:
      #####        KEEP TRACK OF STUFF       #######:
      ##############################################:
      
      keep.rangeU[iter,ttt] = rangeU
      keep.rangeV[iter,ttt] = rangeV
      keep.sigmaU[, iter,ttt] = sigmaU
      keep.sigmaV[, iter,ttt] = sigmaV
      keep.taus1[iter,ttt] = taus1
      keep.taus2[iter,ttt] = taus2
      keep.taue1[iter,ttt] = taue1
      keep.taue2[iter,ttt] = taue2
      keep.A[, iter,ttt] = A
      
      # keep.u1[,,iter]=U1
      # keep.u2[,,iter]=U2
      # keep.v2[,,iter]=V2
      
      keep.Y1.M[,,iter,ttt]=as.matrix(Y1)
      keep.Y2.M[,,iter,ttt]=as.matrix(Y2)
      #print(iter)
    } # End of thin
  }
  print(proc.time()[3] - start)
  
  
  out=list(keep.rangeU,keep.rangeV,keep.sigmaU,keep.sigmaV,keep.taue1,keep.taue2,keep.A,keep.Y1.M,keep.Y2.M)
  names(out)=c('rangeU','rangeV','sigmaU','sigmaV','tau1','tau2','A','Y1.m','Y2.m')
  
}


