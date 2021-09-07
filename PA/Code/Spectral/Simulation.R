library(mvtnorm)
library(MASS)
# Function to compute the MV spatial covariance matrix
# d is the distance matrix between coordinates
# type: t[i] = 1 means AQS; t[i] = 2 mean PA
# theta: parameters
cov.LMC <- function(d,t,theta){
  rho    <- theta[1]
  range1 <- exp(theta[2])
  range2 <- exp(theta[3])
  tau1   <- exp(theta[4])
  tau2   <- exp(theta[5])
  sig1   <- exp(theta[6])
  sig2   <- exp(theta[7])
  S      <- diag(ifelse(t==1,tau1,tau2))
  S[t==1,t==1] <- S[t==1,t==1] + sig1*exp(-d[t==1,t==1]/range1)
  S[t==1,t==2] <- S[t==1,t==2] + sig1*rho*exp(-d[t==1,t==2]/range1)
  S[t==2,t==1] <- S[t==2,t==1] + sig1*rho*exp(-d[t==2,t==1]/range1)
  S[t==2,t==2] <- S[t==2,t==2] + sig2*exp(-d[t==2,t==2]/range2) +   
    sig1*rho*rho*exp(-d[t==2,t==2]/range2)
  
  Sinv   <- solve(S)
  logdet <- determinant(S)$modulus[1]
  out    <- list(S=S,Sinv=Sinv,logdet=logdet)
  return(out)
}

# Get covariance matrix for U1, U2, and V2
# U1 ~ normal(0, sig1 * Sigma11)
# U2 ~ normal(0, sig1 * Sigma22)
# V2 ~ normal(0, sig2 * Sigma22)
# Cov(U1, U2) = sig1 * Sigma12
cov <- function(d,t,theta) {
  rho    <- theta[1]
  range1 <- exp(theta[2])
  range2 <- exp(theta[3])
  tau1   <- exp(theta[4])
  tau2   <- exp(theta[5])
  sig1   <- exp(theta[6])
  sig2   <- exp(theta[7])
  S <- diag(0, length(t))
  S[t==1,t==1] <- sig1*exp(-d[t==1,t==1]/range1)
  S[t==1,t==2] <- sig1*exp(-d[t==1,t==2]/range1)
  S[t==2,t==1] <- sig1*exp(-d[t==2,t==1]/range1)
  S[t==2,t==2] <- sig2*exp(-d[t==2,t==2]/range2)
  return(S)
}


set.seed(123)
n      <- 250
type   <- rbinom(n,1,0.7)+1
s      <- cbind(runif(n),runif(n))
rho    <- 1
range1 <- 0.1
range2 <- 0.3
tau1   <- 0.1
tau2   <- 0.3
sig1   <- 1
sig2   <- 4
theta  <- c(rho,log(range1),log(range2),
            log(tau1),log(tau2),log(sig1),log(sig2))
Sig    <- cov.LMC(as.matrix(dist(s)),type,theta)
Y      <- as.vector(t(chol(Sig$S))%*%rnorm(n))

Y1 <- Y[type == 1]
Y2 <- Y[type == 2]
n1 <- length(Y1)
n2 <- length(Y2)

# Covariance matrix:
S <- cov(as.matrix(dist(s)),type,theta)
S11 <- S[type == 1, type == 1]
S12 <- S[type == 1, type == 2]
S21 <- S[type == 2, type == 1]
S22 <- S[type == 2, type == 2]

# Simulation for U1: 
set.seed(123)
U2 <- as.vector(t(chol(sig1 * S22)) %*% rnorm(n2))
S1 <- S11 - S12 %*% solve(S22) %*% S21
A1 <- S12 %*% solve(S22)
S1inv <- solve(S1)

SigmaU1 <- solve(1/tau1 * diag(1, n1) + 1/sig1 * S1inv)
meanU1 <- SigmaU1 %*% (1/tau1 * Y1 + 1/sig1 * S1inv %*% A1 %*% U2)

U1_sim <- as.vector(t(chol(SigmaU1)) %*% rnorm(n1)) + meanU1
hist(U1_sim)

# Simulation for U2:
# A_l = rho in this setting
set.seed(123)
U1 <- as.vector(t(chol(sig1 * S11)) %*% rnorm(n1))
V2 <- as.vector(t(chol(sig2 * S22)) %*% rnorm(n2))
S2 <- S22 - S21 %*% solve(S11) %*% S12
A2 <- S21 %*% solve(S11)
S2inv <- solve(S2)

SigmaU2 <- solve(1/tau2 * diag(1, n2) + 1/sig1 * S2inv)
meanU2 <- SigmaU2 %*% (1/tau2 * rho * Y2 - rho * V2 + 1/sig1 * S2inv %*% A2 %*% U1)

U2_sim <- as.vector(t(chol(SigmaU2)) %*% rnorm(n2)) + meanU2
hist(U2_sim)

# Simulation for V2:
set.seed(123)
SigmaV2 <- solve(solve(S22) + 1/tau2 * diag(1, n2))
meanV2 <- SigmaV2 %*% (1/tau2 * Y2)

V2_sim <- as.vector(t(chol(SigmaV2)) %*% rnorm(n2)) + meanV2
hist(V2_sim)





















