rm(list = ls())
library(mvtnorm)
library(MASS)
library(truncnorm)
library(invgamma)
########################################################################
##########      Functions
########################################################################
# Get spatial covariance function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}
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
  S[t==2,t==2] <- sig1*exp(-d[t==2,t==2]/range1)
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
n      <- 200
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
tau1   <- 0.1
tau2   <- 0.3
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