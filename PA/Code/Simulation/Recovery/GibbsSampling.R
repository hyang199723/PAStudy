rm(list = ls())
library(mvtnorm)
library(MASS)
library(truncnorm)
library(invgamma)
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Spectral/")
source('SimData.R')

########################################################################
##########      Generate Initial Values for U1, U2, V2, Al, and sigmas
########################################################################
U1_init <- as.vector(rnorm(n1))
U2_init <- as.vector(rnorm(n2))
V2_init <- as.vector(rnorm(n2))
Al_init <- 0
sig1_init <- 2
sig2_init <- 2
tau1_init = 0.1
tau2_init = 0.3
########################################################################
##########   Gibbs Sampler
########################################################################
iter = 2000
## Parameters required to generate samples
S <- covIndividual(dist, type, theta)
S11 <- S$S11; S12 <- S$S12; S21 <- S$S21; S22 <- S$S22
# U1
S1 <- S11 - S12 %*% solve(S22) %*% S21
S1inv <- solve(S1)
A1 <- S12 %*% solve(S22)

U1_sim_all <- matrix(data = 0, nrow = n1, ncol = iter)
# U2
S2 <- S22 - S21 %*% solve(S11) %*% S12
A2 <- S21 %*% solve(S11)
S2inv <- solve(S2)
U2_sim_all <- matrix(data = 0, nrow = n2, ncol = iter)
# V2
V2_sim_all <- matrix(data = 0, nrow = n2, ncol = iter)
# Al
Al_sim_all <- rep(0,iter)
# sig1
sig1_sim_all <- rep(0, iter)
# sig2
sig2_sim_all <- rep(0, iter)
# tau1
tau1_sim_all <- rep(0, iter)
# tau2
tau2_sim_all <- rep(0, iter)

## Sampling
tau1_sim = tau1_init
tau2_sim = tau2_init
U2_sim = U2_init
V2_sim = V2_init
Al_sim = Al_init
sig1_sim = sig1_init
sig2_sim = sig2_init

for (i in 1:iter) {
  # Sample U1
  sigmaU1 <- solve(1/tau1_sim * diag(1, n1) + 1/sig1_sim * S1inv)
  meanU1 <- sigmaU1 %*% (1/tau1_sim * Y1 + 1/sig1_sim * S1inv %*% A1 %*% U2_sim)
  U1_sim <- as.vector(t(chol(sigmaU1)) %*% rnorm(n1)) + meanU1
  U1_sim_all[,i] = U1_sim
  
  # Sample U2
  sigmaU2 <- solve(Al_sim^2/tau2_sim * diag(1, n2) + 1/sig1_sim * S2inv)
  meanU2 <- sigmaU2 %*% (1/tau2_sim * Al_sim * Y2 - 1/tau2_sim * Al_sim * V2_sim + 
                           1/sig1_sim * S2inv %*% A2 %*% U1_sim)
  U2_sim <- as.vector(t(chol(sigmaU2)) %*% rnorm(n2)) + meanU2
  U2_sim_all[,i] = U2_sim
  
  # Sample V2
  sigmaV2 <- solve(1/sig2_sim *solve(Sv1) + 1/tau2_sim * diag(1, n2))
  meanV2 <- sigmaV2 %*% (1/tau2_sim * Y2 - 1/tau2_sim * Al_sim * U2_sim)
  V2_sim <- as.vector(t(chol(sigmaV2)) %*% rnorm(n2)) + meanV2
  V2_sim_all[,i] = V2_sim
  
  # Sample Al
  sigmaAl <- solve(1/tau2_sim * t(U2_sim) %*% U2_sim + 1/5)
  meanAl <- sigmaAl %*% (1/tau2_sim * t(Y2) %*% U2_sim - 1/tau2_sim * t(V2_sim) %*% U2_sim + 2/5)
  Al_sim <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl, sd=sqrt(sigmaAl))
  Al_sim_all[i] = Al_sim
  
  # Sample sig1
  U_sim <- as.vector(append(U1_sim, U2_sim))
  a <- (n2+n1)/2 + 1
  b <- (t(U_sim) %*% solve(S$S) %*% U_sim) / 2 + 1
  sig1_sim <- rinvgamma(1, a, b)
  sig1_sim_all[i] = sig1_sim
  
  # Sample sig2
  a <- n2/2 + 1
  b <- (t(V2_sim) %*% solve(Sv1) %*% V2_sim) / 2 + 1
  sig2_sim <- rinvgamma(1, a, b)
  sig2_sim_all[i] = sig2_sim
  
  # Sample tau1
  a <- n1/2 + 1
  b <- t(Y1 - U1) %*% (Y1 - U1) / 2 + 1
  tau1_sim <- rinvgamma(1, a, b)
  tau1_sim_all[i] = tau1_sim
  
  # Sample tau2
  a <- n2/2 + 1
  b <- t(Y2 - Al*U2 - V2) %*% (Y2 - Al*U2 - V2) / 2 + 1
  tau2_sim <- rinvgamma(1, a, b)
  tau2_sim_all[i] = tau2_sim
}
### Plot results
par(mfrow=c(3,3))
index <- c(2,20,30)
for (j in index) {
  plot(U1_sim_all[j,],type="l",
       xlab="MCMC iteration",ylab="Sample",
       main=paste('U1 station',j))
  abline(U1[j],0,col=2,lwd=2) 
}
for (j in index) {
  plot(U2_sim_all[j,],type="l",
       xlab="MCMC iteration",ylab="Sample",
       main=paste('U2 station',j))
  abline(U2[j],0,col=2,lwd=2) 
}
for (j in index) {
  plot(V2_sim_all[j,],type="l",
       xlab="MCMC iteration",ylab="Sample",
       main=paste('V2 station',j))
  abline(V2[j],0,col=2,lwd=2) 
}
plot(sig1_sim_all,type="l",
     xlab="MCMC iteration",ylab="Sample",
     main='sigma1')
abline(sig1,0,col=2,lwd=2)
plot(sig2_sim_all,type="l",
     xlab="MCMC iteration",ylab="Sample",
     main='sigma2')
abline(sig2,0,col=2,lwd=2)
plot(tau1_sim_all,type="l",
     xlab="MCMC iteration",ylab="Sample",
     main='tau1')
abline(tau1,0,col=2,lwd=2)
plot(tau2_sim_all,type="l",
     xlab="MCMC iteration",ylab="Sample",
     main='tau2')
abline(tau2,0,col=2,lwd=2)






