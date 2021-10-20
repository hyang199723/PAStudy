rm(list=ls())
setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Spectral/Hongjian/ExtraFunctions.R')
library(fields) 
library(geoR)
library(truncnorm)
library(tidyverse)
library(mvtnorm)


# Put test set as NAs 
# do cross-validation leave-one out cross validation 
# Bias variance squared error 95% converage
# Average prediction interval
# Compare it with or without PA data
# Plot correlation v.s. frequency


# 

# Plot Al vs frequency
# plot variance vs frequency

# compare this algorithm with normal kriging


# Read in data
PA_data <- read.csv("Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv")
FRM_data <- read.csv("Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv")
# Convert timestamp
PA_data$Timestamp <- as.POSIXct(PA_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
FRM_data$Timestamp <- as.POSIXct(FRM_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
# Start with one week data to check performance

# No PA 
# 
start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-03 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7



pa <- subset(PA_data, (Timestamp >= start) & (Timestamp <= end))
frm <- subset(FRM_data, (Timestamp >= start) & (Timestamp <= end))
# Get data to desired format
paTS <- pivot_wider(pa, names_from = Timestamp, values_from = PM25)
frmTS <- pivot_wider(frm, names_from = Timestamp, values_from = PM25)
# Record locations of PA and FRM stations
s1 <- frmTS[, 1:2]
s2 <- paTS[, 1:2]
paTS <- paTS[, -c(1:2)]
frmTS <- frmTS[, -c(1:2)]
Y1 = data.frame(frmTS)
Y2 = data.frame(paTS)

# set initial values 
mean_range=0
sd_range=1
mean_var=0
sd_var=1
mean_rho=0
sd_rho=10
iters=3000
burn=1000
thin=1
update=10


# Bookkeeping
n1       <- nrow(Y1)
n2       <- nrow(Y2)
nt       <- ncol(Y1)
m1       <- is.na(Y1)
m2       <- is.na(Y2)
d        <- as.matrix(dist(rbind(s1,s2)))
dv2=as.matrix(dist(s2))
const    <- nt/2# ?

# Initial values
beta1 <- mean(colMeans(Y1,na.rm=TRUE))
beta2 <- mean(colMeans(Y2,na.rm=TRUE))
Y1[m1] <- beta1
Y2[m2] <- beta2
U1     <- matrix(0,n1,nt)
U2     <- matrix(0,n2,nt)
V2     <- matrix(0,n2,nt)
rangeU  <- exp(mean_rho)
rangeV   <- exp(mean_rho)
lrangeU  <- log(rangeU)
lrangeV   <-log(rangeV)

taue1  <- mean(sapply(Y1, sd)^2)
taue2  <- mean(sapply(Y2, sd)^2)
sigmaU   <- rep(taue1,nt)
sigmaV   <- rep(taue1,nt)
A      <- rep(1,nt)
Z1     <- matrix(0,n1,nt)
Z2     <- matrix(0,n2,nt)
Ys1    <- matrix(0,n1,nt)
Ys2    <- matrix(0,n2,nt)

# Keep track of stuff

keep_theta  <- array(0,dim=c(iters,9,nt))
keep.u1= array(NA,dim=c(iters,n1,nt))
keep.u2= array(NA,dim=c(iters,n2,nt))
keep.v2= array(NA,dim=c(iters,n2,nt))
keep.Y1.M= array(NA,dim=c(iters,n1,nt))
keep.Y2.M= array(NA,dim=c(iters,n2,nt))

start = proc.time()[3]

# 
for(iter in 1:iters){
  
  for(ttt in 1:thin){
    
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
    taue1 <- 1/rgamma(1,n1*nt/2+0.01,sum((Y1-beta1-Z1)^2)/2+0.01)
    taue2 <- 1/rgamma(1,n2*nt/2+0.01,sum((Y2-beta2-Z2)^2)/2+0.01)
    # 
    ##############################################:
    ####     Transform to spectral land      #####:
    ##############################################:
    
    for(i in 1:n1){
      Ys1[i,] <- fft_real(as.numeric(Y1[i,]-beta1))
    }
    for(i in 1:n2){
      Ys2[i,] <- fft_real(as.numeric(Y2[i,]-beta2))
    }
    taus1 <- const*taue1
    taus2 <- const*taue2
    
    ##############################################:
    ####      LMC TERMS (spectral space)     #####:
    ##############################################:
    # We first keep range the same for all spectrum
    
    eigU = invU(d,n1,n2,rangeU)
    S11 = eigU$S11; S12 <- eigU$S12; S21 <- eigU$S21; S22 <- eigU$S22
    S1 = S11 - eigU$A12 %*% S21
    ES1=eigen(S1)
    S1_G   = ES1$vectors
    S1_D   = ES1$values
    S1inv = S1_G%*%diag(1/S1_D)%*%t(S1_G)
    A1 = S12 %*% solve(S22)
    # S1inv = G * diag(1/D) * G'
    # inverse = G * [1/taus1 + 1/sigmaU * diag(1/D)]^{-1} * G'
    
    S11inv=solve(S11)
    S2 = S22 - S21 %*% S11inv %*% S12
    A2 = S21 %*% S11inv
    ES2 = eigen(S2)
    S2_G = ES2$vectors
    S2_D = ES2$values
    S2inv = S2_G%*%diag(1/S2_D)%*%t(S2_G)
    
    # G: eigenvectors
    # D: eigenvalues
    # Q: inverse matrix
    eigV  = invV(dv2,rangeV)
    V_G = eigV$G
    V_D = eigV$D
    
    for (r in 1:nt) # for each spectral 
    {
      # # Sample U1
      newDiag = diag(1/taus1 * diag(1, n1) + 1/sigmaU[r] * diag(1/S1_D))
      sigmaU1 <- S1_G %*% diag(1/newDiag) %*% t(S1_G)
      # S1_G
      meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1inv %*% A1 %*% U2[,r])
      U1[,r] <- as.vector(t(chol(sigmaU1)) %*% rnorm(n1)) + meanU1 # Generate U matrix with S1_G %*% Normal
      
      # # Sample U2
      newDiag = diag(A[r]^2/taus2 * diag(1, n2) + 1/sigmaU[r] * diag(1/S2_D))
      sigmaU2 <- S2_G %*% diag(1/newDiag) %*% t(S2_G)

      meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 * A[r] * V2[,r] +
                               1/sigmaU[r] * S2inv %*% A2 %*% U1[,r])
      U2[,r] <- as.vector(t(chol(sigmaU2)) %*% rnorm(n2)) + meanU2
      
      # # Sample V2
      newDiag = diag(1/sigmaV[r] * diag(1/V_D) + 1/taus2 * diag(1, n2))
      sigmaV2 = V_G %*% diag(1/newDiag) %*% t(V_G)

      meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
      V2[,r] <- as.vector(t(chol(sigmaV2)) %*% rnorm(n2)) + meanV2
      
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
    # Ru should be a matrix and use data from all spectrum
    # Sweep operation to get rid of variance
    
    Ru = as.vector(U_sim)/sqrt(sigmaU[r])
    Rv = as.vector(V2[,r])/sqrt(sigmaV[r])
    
    # range1
    Ms=exp_corr(d,range=exp(lrangeU))
    curll = dmvnorm(Ru,rep(0,n1+n2),Ms,log=TRUE)
    canrange1 = rnorm(1,lrangeU,0.5)
    canM = exp_corr(d,range=exp(canrange1))
    canll = dmvnorm(Ru,rep(0,n1+n2),canM,log=TRUE)
    
    MH1 <- canll-curll+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
    
    if (log(runif(1))<MH1)
    {
      lrangeU=canrange1
    }
    
    # range2
    Ss=exp_corr(dv2,range = exp(lrangeV))
    curll2 = dmvnorm(Rv,rep(0,n2),Ss,log=TRUE)
    canrange2 = rnorm(1,lrangeV,0.5)
    canS = exp_corr(dv2,range=exp(canrange2))
    canll2 = dmvnorm(Rv,rep(0,n2),canS,log=TRUE)
    
    MH2 <- canll2-curll2+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
    
    if (log(runif(1))<MH2)
    {
      lrangeV=canrange2
    }
    rangeU = exp(lrangeU)
    rangeV = exp(lrangeV)
    
  } # end thin
  
  ##############################################:
  #####        KEEP TRACK OF STUFF       #######:
  ##############################################:
  
  
  keep_theta[iter, 1,]  <- rep(rangeU,nt)
  keep_theta[iter, 2,]  <- rep(rangeV,nt)
  keep_theta[iter, 3,]  <- sigmaU
  keep_theta[iter, 4,]  <- sigmaV
  keep_theta[iter, 5,]  <- rep(taue1,nt) # because we have the same tau for all spectrums 
  keep_theta[iter, 6,]  <- rep(taue2,nt)
  keep_theta[iter, 7,]  <- A
  keep_theta[iter, 8,]  <- rep(beta1,nt)
  keep_theta[iter, 9,]  <- rep(beta2,nt)
  
  keep.u1[iter,,]=U1
  keep.u2[iter,,]=U2
  keep.v2[iter,,]=V2
  
  keep.Y1.M[iter,,]=as.matrix(Y1)
  keep.Y2.M[iter,,]=as.matrix(Y2)
}
proc.time()[3] - start

### Results 

# create data frame 
P=c("rangeu","rangev","sigmaU","sigmaV","taue1","taue2","A",'beta1','beta2')
rv=list(rep(rangeU,nt),rep(rangeV,nt),sigmaU,sigmaV,rep(taue1,nt),rep(taue2,nt),A,
        rep(beta1,nt),rep(beta2,nt))
values=c()
param=c()
n.iter=c()
time=c()
realvalues=c()
for(i in 1:9)
{
  for (j in 1:nt)
  {
    values=c(values,keep_theta[,i,j])
    param=c(param,rep(P[i],iters))
    n.iter=c(n.iter,rep(1:iters))
    time=c(time,rep(j,iters))
    realvalues=c(realvalues,rep(rv[[i]][j],iters))
    
  }
}
# 19 timestamps - 10 min
# 43 timestamps - 19
# 67 TS - 28 min


mcmc.results=data.frame(n.iter,values,param,time,realvalues)
mcmc.results$n.iter=as.numeric(mcmc.results$n.iter)

#al 
ggplot(mcmc.results %>% filter(param=='A',time%in%seq(1:10)))+geom_line(aes(x=n.iter,y=values))+
  facet_wrap(~time, scales = "free")+theme_bw()+
  ggtitle('Al values')
#tau1
ggplot(mcmc.results %>% filter(param=='taue1'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+ggtitle('tau1 values')
#tau2
ggplot(mcmc.results %>% filter(param=='taue2'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+ggtitle('tau2 values')

#sigmau
ggplot(mcmc.results %>% filter(param=='sigmaU'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+ggtitle('sigma1 values')+
  facet_wrap(~time,scales = "free")

#sigmav
ggplot(mcmc.results %>% filter(param=='sigmaV'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+ggtitle('sigma2 values')+
  facet_wrap(~time,scales = "free")

#phou
ggplot(mcmc.results %>% filter(param=='rangeu'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+ggtitle('phou values')

ggplot(mcmc.results %>% filter(param=='rangev'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+
  ggtitle('phov')


