rm(list = ls())
library(mvtnorm)
library(MASS)
library(truncnorm)
library(invgamma)
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Spectral/")
source('SimRecovery/SimData.R')

################################################
########## Functions
################################################
# Spatial covariance function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}
# Log likelihood function of posterior distribution
# Y: observed values (multi-normal)
# n: Number of observations
# range: range1 or range2, current values (normal)
# distance: distance matrix
log_post <- function(Y, n, range, distance) {
  covariance = exp_corr(distance, range)
  like = dmvnorm(Y, rep(0, n), covariance, log=TRUE)
  prior = dnorm(range, log=TRUE)
  return(like + prior)
}


################################################
########## MH within Gibbs
################################################
iters = 10000
# Generate Initial Values for U1, U2, V2, Al, and sigmas
U1_init <- as.vector(rnorm(n1))
U2_init <- as.vector(rnorm(n2))
V2_init <- as.vector(rnorm(n2))
Al_init <- 0.2
sig1_init <- 2
sig2_init <- 2
tau1_init = 0.1
tau2_init = 0.3
########################################################################
##########   Gibbs Sampler
########################################################################
# U1
U1_sim_all <- matrix(data = 0, nrow = n1, ncol = iters)
# U2
U2_sim_all <- matrix(data = 0, nrow = n2, ncol = iters)
# V2
V2_sim_all <- matrix(data = 0, nrow = n2, ncol = iters)
# Al
Al_sim_all <- rep(0,iters)
# sig1
sig1_sim_all <- rep(0, iters)
# sig2
sig2_sim_all <- rep(0, iters)
# tau1
tau1_sim_all <- rep(0, iters)
# tau2
tau2_sim_all <- rep(0, iters)

## Sampling
tau1_sim = tau1_init
tau2_sim = tau2_init
U1_sim = U1_init
U2_sim = U2_init
V2_sim = V2_init
Al_sim = Al_init
sig1_sim = sig1_init
sig2_sim = sig2_init
U_sim = c(U1_sim, U2_sim)

### Metropolis Hasting 
# Bookkeeping
keep.range1 = rep(NA, iters)
keep.range2 = rep(NA, iters)

# Initial values
#Ru = as.vector(u)/sqrt(sigmau)
#Rv = as.vector(v2)/sqrt(sigmav)
#U = U / sqrt(sig1)
#V2 = V2 / sqrt(sig2) #?

# Get distance matrix:
dist_full = dist
dist22 = dist[(n1+1):(n1 + n2), (n1+1):(n1 + n2)]

currange1 = 0.5
currange2 = 0.5
for(i in 1:iters){
  # range1
  canrange1 = rnorm(1, currange1, 0.5)
  currU = U_sim/sqrt(sig1_sim)
  logR1 <- log_post(currU, n1 + n2, canrange1, dist_full) - log_post(currU, n1 + n2, currange1, dist_full) 
  
  if (log(runif(1)) < logR1) {
    currange1 = canrange1
  }
  keep.range1[i]  <- currange1
  
  # range2
  canrange2 = rnorm(1, currange2, 0.5)
  currV = as.vector(V2_sim/sqrt(sig2_sim))
  logR2 <- log_post(currV, n2, canrange2, dist22) - log_post(currV, n2, currange2, dist22) 
  
  if (log(runif(1)) < logR2)
  {
    currange2 = canrange2
  }
  keep.range2[i]  <- currange2
  
  ## Parameters required to generate Gibbs samples
  currS <- exp_corr(dist_full, currange1)
  S11 <- currS[1:n1, 1:n1]
  S12 <- currS[1:n1, (n1+1):(n1+n2)]
  S21 <- currS[(n1+1):(n1+n2), 1:n1]
  S22 <- currS[(n1+1):(n1+n2), (n1+1):(n1+n2)]
  # U1 and its eigens
  S1 <- S11 - S12 %*% solve(S22) %*% S21
  eigS1 <- eigen(S1)
  S1_G = eigS1$vectors
  S1_D = eigS1$values
  S1inv <- S1_G %*% diag(1/S1_D) %*% t(S1_G)
  
  A1 <- S12 %*% solve(S22)
  # U2
  S2 <- S22 - S21 %*% solve(S11) %*% S12
  A2 <- S21 %*% solve(S11)
  eigS2 = eigen(S2)
  S2_G = eigS2$vectors
  S2_D = eigS2$values
  
  S2inv <- S2_G %*% diag(1/S2_D) %*% t(S2_G)
  
  # Sample U1
  newDiag = 1/tau1_sim + 1/sig1_sim * 1/S1_D
  sigmaU1 <- S1_G %*% diag(1/newDiag) %*% t(S1_G)
  meanU1 <- sigmaU1 %*% (1/tau1_sim * Y1 + 1/sig1_sim * S1inv %*% A1 %*% U2_sim)
  U1_sim <- meanU1+S1_G%*%rnorm(n1, 0, 1/sqrt(newDiag))
  U1_sim_all[,i] = U1_sim
  
  # Sample U2
  newDiag = Al_sim^2/tau2_sim + 1/sig1_sim * 1/S2_D
  sigmaU2 <- S2_G %*% diag(1/newDiag) %*% t(S2_G)
  sigmaU2 <- solve(Al_sim^2/tau2_sim * diag(1, n2) + 1/sig1_sim * S2inv)
  meanU2 <- sigmaU2 %*% (1/tau2_sim * Al_sim * Y2 - Al_sim * V2_sim + 
                           1/sig1_sim * S2inv %*% A2 %*% U1_sim)
  U2_sim <- meanU2+S2_G%*%rnorm(n2,0,1/sqrt(newDiag))
  U2_sim_all[,i] = U2_sim
  
  # Sample V2
  currSv <- exp_corr(dist22, currange2)
  sigmaV2 <- solve(1/sig2_sim * solve(currSv) + 1/tau2_sim * diag(1, n2))
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
  b <- (t(U_sim) %*% solve(currS) %*% U_sim) / 2 + 1
  sig1_sim <- rinvgamma(1, a, b)
  sig1_sim_all[i] = sig1_sim
  
  # Sample sig2
  a <- n2/2 + 1
  b <- (t(V2_sim) %*% solve(currSv) %*% V2_sim) / 2 + 1
  sig2_sim <- rinvgamma(1, a, b)
  sig2_sim_all[i] = sig2_sim
  
  # Sample tau1
  a <- n1/2 + 1
  b <- t(Y1 - U1_sim) %*% (Y1 - U1_sim) / 2 + 1
  tau1_sim <- rinvgamma(1, a, b)
  tau1_sim_all[i] = tau1_sim
  
  # Sample tau2
  a <- n2/2 + 1
  b <- t(Y2 - Al_sim*U2_sim - V2_sim) %*% (Y2 - Al_sim*U2_sim - V2_sim) / 2 + 1
  tau2_sim <- rinvgamma(1, a, b)
  tau2_sim_all[i] = tau2_sim
}
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


plot(x = 1:iters, y = keep.range1, 'l', main = 'range1')
abline(h = 0.1, col = 'red')

plot(x = 1:iters, y = keep.range2, 'l', main = 'range2')
abline(h = 0.3, col = 'red')
