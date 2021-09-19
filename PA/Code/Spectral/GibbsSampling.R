rm(list = ls())
library(mvtnorm)
library(MASS)
library(truncnorm)
library(invgamma)
########################################################################
##########      Functions
########################################################################
# Get covariance matrix for U1, U2
# d: distance matrix between coords
# t: type, either 1 or 2
# theta: parameters
covU <- function(d,t,theta) {
  range1 <- exp(theta[2])
  sig1   <- exp(theta[6])
  S <- diag(0, length(t))
  S[t==1,t==1] <- sig1*exp(-d[t==1,t==1]/range1)
  S[t==1,t==2] <- sig1*exp(-d[t==1,t==2]/range1)
  S[t==2,t==1] <- sig1*exp(-d[t==2,t==1]/range1)
  S[t==2,t==2] <- sig1*exp(-d[t==2,t==2]/range1) # Should it be range1 or range2?
  return(S)
}

# Get covariance matrix for V2, Sigma_{v22}
# d: distance matrix between coords
# t: type, either 1 or 2
# theta: parameters
covV <- function(d,t,theta) {
  # We only need type = 2
  d <- d[t == 2, t== 2]
  t <- t[t == 2]
  range2 <- exp(theta[3])
  sig2   <- exp(theta[7])
  S <- diag(0, length(t))
  S <- sig2*exp(-d/range2)
  return(S)
}

# Function to compute Sigma11, Sigma12, Sigma21, and Sigma22
# d: distance matrix between coords
# t: type, either 1 or 2
# theta: parameters
# # theta  <- c(rho,log(range1),log(range2),log(tau1),log(tau2),log(sig1),log(sig2))
covIndividual <- function(d, t, theta) {
  range1 <- exp(theta[2])
  S <- diag(0, length(t))
  S[t==1,t==1] <- exp(-d[t==1,t==1]/range1)
  S[t==1,t==2] <- exp(-d[t==1,t==2]/range1)
  S[t==2,t==1] <- exp(-d[t==2,t==1]/range1)
  S[t==2,t==2] <- exp(-d[t==2,t==2]/range1)
  S11 <- S[t==1,t==1]
  S12 <- S[t==1,t==2]
  S21 <- S[t==2,t==1]
  S22 <- S[t==2,t==2]
  return(list(S=S,S11=S11, S12=S12, S21=S21, S22=S22))
}

########################################################################
##########      Simulate data
########################################################################
set.seed(123)
# Total number of observations
n      <- 250
# 70% PA(Type 2), 30% AQS (Type 1)
type   <- rbinom(n,1,0.7)+1
# Ordered, easier for simulation
type <- type[order(type)]
# Generate coords
s      <- cbind(runif(n),runif(n))
# Distance matrix
dist <- as.matrix(dist(s))
# Parameters
# \sigma_{ul} = \sig1, #\sigma_{vl} = \sig2
rho    <- 1 # rho was not used, but left here for legacy reason
range1 <- 0.1
range2 <- 0.3
sig1 <- 1
sig2 <- 2
theta  <- c(rho,log(range1),log(range2),
            log(tau1),log(tau2),log(sig1),log(sig2))
n1 <- sum(type == 1)
n2 <- sum(type == 2)

# Covariance matrix for U1, U2
Su <- covU(dist,type,theta)
Su1 <- covIndividual(dist,type,theta)
# Covariance matrix for V2
Sv <- covV(dist,type,theta)
Sv1 <- Sv/sig2


# Generate U1, U2, V2, Al
# Generate U1, U2
# The first n1 values are U1, remaining are U2
set.seed(129)
U <- as.vector(t(chol(Su1$S)) %*% rnorm(n1 + n2, mean = 0, sd = sqrt(sig1)))
U1 <- U[1:n1]
U2 <- U[n1+1:n2]
if (length(U1) != n1 | length(U2) != n2) {stop('length do not match')}

# Generate V2
set.seed(130)
V2 <- as.vector(t(chol(Sv1)) %*% rnorm(n2, mean = 0, sd = sqrt(sig2)))
# Generate A_{l}
Al <- 0.8

# Generate Y1 and Y2
set.seed(131)
Y1 <- U1 + as.vector(rnorm(n1, mean = 0, sd = sqrt(tau1)))
Y2 <- Al * U2 + V2 + as.vector(rnorm(n2, mean = 0, sd = sqrt(tau2)))

########################################################################
##########      Generate Initial Values for U1, U2, V2, Al, and sigmas
########################################################################
set.seed(30)
U1_init <- as.vector(rnorm(n1))
set.seed(31)
U2_init <- as.vector(rnorm(n2))
set.seed(32)
V2_init <- as.vector(rnorm(n2))
Al_init <- 0
sig1_init <- 2
sig2_init <- 2

########################################################################
##########   Gibbs Sampler
########################################################################
iter <- 3000
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
##### Note: tau1 and tau2 also needs updating
tau1_sim = tau1
tau2_sim = tau2
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
