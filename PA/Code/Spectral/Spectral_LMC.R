source('ExtraFunctions.R')
source('simAllTS.R')
#######################################################
#
# INPUTS
#
# Y1         := n1 x nt-vector of observations
# Y2         := n2 x nt-vector of observations
# s1         := n1x2 matrix of spatial coordinates
# s2         := n2x2 matrix of spatial coordinates
# mean_range := log(spatial range) ~ N(mean_range,sd_range)
# sd_range
# mean_var   := log(variance parameters) ~ N(mean_var,sd_var)
# sd_var
# mean_rho   := rho ~ N(mean_rho,sd_rho)
# sd_rho
#
# iters      := number of MCMC iterations
# burn       := length of burn-in
#
# Model
#
# Y1[i,t]    = beta1 + U1[i,t] + e1[i,t]
# Y2[i,t]    = beta2 + A*U1[i,t] + U2[i,t] + e2[i,t]
# Var(Z) = sigU
# Var(Z) = sigV
# Var(e1) = tau1
# Var(e2) = tau2
# Cor(Z1(s1),Z1(s2)) = exp(-d12/range1)
# Cor(Z2(s1),Z2(s2)) = exp(-d12/range2)
#
#######################################################
library(fields) 
library(geoR)
library(truncnorm)

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

# The main MCMC function
# LMC <- function(Y1,Y2,s1,s2,
#                 mean_range=0,sd_range=1,
#                 mean_var=0,sd_var=1,
#                 mean_rho=0,sd_rho=10,
#                 iters=5000,burn=1000,thin=1,update=10){
  
# Bookkeeping
n1       <- nrow(Y1)
n2       <- nrow(Y2)
nt       <- ncol(Y1)
m1       <- is.na(Y1)
m2       <- is.na(Y2)
d        <- as.matrix(dist(rbind(s1,s2)))
const    <- nt# 1/?
  
# Initial values
  
beta1 <- mean(Y1,na.rm=TRUE)
beta2 <- mean(Y2,na.rm=TRUE)
Y1[m1] <- beta1
Y2[m2] <- beta2
U1     <- matrix(0,n1,nt)
U2     <- matrix(0,n2,nt)
V2     <- matrix(0,n2,nt)
rangeU  <- rep(exp(mean_rho),nt)
rangeV   <- rep(exp(mean_rho),nt)
taue1  <- var(as.vector(Y1)) #1/?
taue2  <- var(as.vector(Y2)) #1/?
sigmaU   <- rep(taue1,nt)
sigmaV   <- rep(taue1,nt)
A      <- rep(1,nt)
Z1     <- matrix(0,n1,nt)
Z2     <- matrix(0,n2,nt)
Ys1    <- matrix(0,n1,nt)
Ys2    <- matrix(0,n2,nt)
# eigU   <- invU(d,n1,n2,rhoU)
# eigV   <- invV(d,rhoV)
  
  
# Keep track of stuff
  
keep_theta  <- array(0,dim=c(iters,9,nt))
keep.u1= array(NA,dim=c(iters,n1,nt))
keep.u2= array(NA,dim=c(iters,n2,nt))
keep.v2= array(NA,dim=c(iters,n2,nt))



# set to the real values 
  
lrangeU  <- lrangeu
lrangeV   <-lrangev
taue1  <- tau1 #1/?
taue2  <- tau2 #1/?
sigmaU   <- sigmau
sigmaV   <- sigmav

# GO!!!
pb <- txtProgressBar(min = 0, max = iters, style = 3)
  
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
      
      Y1[m1] <- rnorm(sum(m1),beta1+Z1[m1],sqrt(taue1)) #1/?
      Y2[m2] <- rnorm(sum(m2),beta2+Z2[m2],sqrt(taue2)) #1/?
      
      ##############################################:
      ####      MEAN PARAMETERS (real space)   #####:
      ##############################################:
      
      # full conditional for beta1 and beta2 ...
      #VVV   <- (taue1*n1*nt + 0.01)
      #MMM   <- taue1*sum(Y1-Z1) 
      #beta1 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      beta1=beta.1
      
      # VVV   <- (taue2*n2*nt + 0.01)
      # MMM   <- taue2*sum(Y2-Z2) 
      #beta2 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      beta2=beta.2
      
      ##############################################:
      ####     ERROR VARIANCES (real space)    #####:
      ##############################################:
      # full conditionals for taue1 and taue2 
       taue1 <- 1/rgamma(1,n1*nt/2+0.01,sum((Y1-beta1-Z1)^2)/2+0.01)
       taue2 <- 1/rgamma(1,n2*nt/2+0.01,sum((Y2-beta2-Z2)^2)/2+0.01)
      
      ##############################################:
      ####     Transform to spectral land      #####:
      ##############################################:
      
      for(i in 1:n1){Ys1[i,] <- fft_real(Y1[i,]-beta1)}
      for(i in 1:n2){Ys2[i,] <- fft_real(Y2[i,]-beta2)}
      taus1 <- const*taue1
      taus2 <- const*taue2
      
      ##############################################:
      ####      LMC TERMS (spectral space)     #####:
      ##############################################:
      
      for (r in 1:nt) # for each spectral 
      {
        
        #eigU <-invU(d,n1,n2,rangeU[r]) # is it better way to compute this?
        S=covIndividual(d, type, lrangeU[r])
        S11 <- S$S11; S12 <- S$S12; S21 <- S$S21; S22 <- S$S22
        # # Sample U1
        S1 <- S11 - S12 %*% solve(S22) %*% S21
        S1inv <- solve(S1)
        A1 <- S12 %*% solve(S22)
        sigmaU1 <- solve(1/taus1 * diag(1, n1) + 1/sigmaU[r] * S1inv)
        meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1inv %*% A1 %*% U2[,r])
        U1[,r] <- as.vector(t(chol(sigmaU1)) %*% rnorm(n1)) + meanU1
        
        # # Sample U2
        S2 <- S22 - S21 %*% solve(S11) %*% S12
        A2 <- S21 %*% solve(S11)
        S2inv <- solve(S2)
        
        sigmaU2 <- solve(A[r]^2/taus2 * diag(1, n2) + 1/sigmaU[r] * S2inv)
        meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 * A[r] * V2[,r] +
                                 1/sigmaU[r] * S2inv %*% A2 %*% U1[,r])
        
        U2[,r] <- as.vector(t(chol(sigmaU2)) %*% rnorm(n2)) + meanU2
        
        # # Sample V2
        Sv <- covV2(d,type,lrangeV[r])
        
        sigmaV2 <- solve(1/sigmaV[r] *solve(Sv) + 1/taus2 * diag(1, n2))
        meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
        V2[,r] <- as.vector(t(chol(sigmaV2)) %*% rnorm(n2)) + meanV2
        
        # # Sample Al
        sigmaAl <- solve(1/taus2 * t(U2[,r]) %*% U2[,r] + 1/5)
        meanAl <- sigmaAl %*% (1/taus2 * t(Ys2[,r]) %*% U2[,r] - 1/taus2 * t(V2[,r]) %*% U2[,r] + 0.8/5)
        A[r] <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl, sd=sqrt(sigmaAl))
       
      }
      
      
      ###################################################:
      ####  COVARIANCE PARAMETERS (spectral domain) #####:
      ###################################################:
      
      # tauU ... bV are fixed here but should be updated
    #   tauA <- muA  <- 0
    #   aU   <-  bU <- aV <- bV <- 0.1
    #   for(t in 1:nt){
    #     VVV     <- taus2*sum(U2[,t]^2)+tauA
    #     MMM     <- taus2*sum((Ys2[,t]-V2[,t])*U2[,t])+muA*tauA
    #     A[t]    <- rnorm(1,MMM/VVV,1/sqrt(VVV))   
    #     
    #     U       <- c(U1[,t],U2[,t])
    #     tauU[t] <- rgamma(1,(n1+n2)/2+aU,t(U)%*%E$inv%*%U/2+bU)
    #     tauV[t] <- rgamma(1,n2/2+aV,t(V2[,t])%*%E$inv2%*%V2[,t]/2+bV)
    #   }
    #   
    #   
     } # end thin
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    
    keep_theta[iter, 1,]  <- exp(lrangeU)
    keep_theta[iter, 2,]  <- exp(lrangeV)
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

    ##############################################:
    #####       PLOT RESULTS SO FAR        #######:
    ##############################################:
    
    # if(iter%%update==0){
    #   par(mfrow=c(3,3))
    #   for(j in 1:2){
    #     plot(keep_beta[1:iter,j],type="l",
    #          xlab="MCMC iteration",ylab="Sample",
    #          main=colnames(keep_beta)[j])
    #   }
    #   for(j in 1:q){
    #     plot(keep_theta[1:iter,j],type="l",
    #          xlab="MCMC iteration",ylab="Sample",
    #          main=colnames(keep_theta)[j])
    #   }
    # }
    # 
    
  }   
close(pb) 
  #out <- list(beta=keep_beta,theta=keep_theta,
  #            acc_rate=acc/att,time=tock-tick)
  #return(out)}



#### plots
# create data frame 
P=c("rangeu","rangev","sigmaU","sigmaV","taue1","taue2","A",'beta1','beta2')
rv=list(rangeu,rangev,sigmau,sigmav,rep(tau1,nt),rep(tau2,nt),al,
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



mcmc.results=data.frame(n.iter,values,param,time,realvalues)
mcmc.results$n.iter=as.numeric(mcmc.results$n.iter)


ggplot(mcmc.results %>% filter(param=='A',time%in%seq(1:10)))+geom_line(aes(x=n.iter,y=values))+
  facet_wrap(~time)+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))


ggplot(mcmc.results %>% filter(param=='taue1',time%in%seq(1:10)))+geom_line(aes(x=n.iter,y=values))+
  facet_wrap(~time)+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))

ggplot(mcmc.results %>% filter(param=='taue2',time%in%seq(1:10)))+geom_line(aes(x=n.iter,y=values))+
  facet_wrap(~time)+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))



par(mfrow=c(3,4),mar=c(1.3,1.8,1.3,1.3))
for (i in 1:10)
{
  plot(keep_theta[,7,i],type='l')
  abline(h=al[i],col='red')
  
}



s=80
tt=2
plot(keep.u1[,s,tt],type='l')
abline(h=U1[s,tt],col='red')


s=58
tt=2
plot(keep.u2[,s,tt],type='l')
abline(h=U2[s,tt],col='red')


s=25
tt=3
plot(keep.v2[,s,tt],type='l')
abline(h=V2[s,tt],col='red')
