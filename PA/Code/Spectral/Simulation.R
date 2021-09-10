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
  S[t==2,t==2] <- sig1*exp(-d[t==2,t==2]/range1)
  return(S)
}

# Get covariance matrix for V2
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
# Parameters
# \sigma_{ul} = \sig1, #\sigma_{vl} = \sig2
rho    <- 1 # rho was not used, but left here for legacy reason
range1 <- 0.1
range2 <- 0.3
tau1   <- 0.1
tau2   <- 0.3
set.seed(125); sig1 <- rinvgamma(1, 1, 1) # ~invgamma(1,1), sig1 = 0.726
set.seed(124); sig2 <- rinvgamma(1, 1, 0.5) # ~invgamma(1,0.5), sig2 = 1.678
theta  <- c(rho,log(range1),log(range2),
            log(tau1),log(tau2),log(sig1),log(sig2))
n1 <- sum(type == 1)
n2 <- sum(type == 2)

# Covariance matrix for U1, U2
Su <- covU(as.matrix(dist(s)),type,theta)
# Covariance matrix for V2
Sv <- covV(as.matrix(dist(s)),type,theta)

# Generate U1, U2, V2, Al
# Generate U1, U2
# The first n1 values are U1, remaining are U2
set.seed(123)
U <- as.vector(t(chol(Su)) %*% rnorm(n1 + n2))
U1 <- U[1:n1]
U2 <- U[n1+1:n2]
if (length(U1) != n1 | length(U2) != n2) {stop('length do not match')}
# Generate V2
set.seed(123)
V2 <- as.vector(t(chol(Sv)) %*% rnorm(n2))
# Generate A_{l}
set.seed(123); Al <- rtruncnorm(1, a=0, b=+Inf, mean=2, sd=2) # 0.879

# Generate Y1 and Y2
set.seed(123)
Y1 <- U1 + as.vector(rnorm(n1, mean = 0, sd = sqrt(tau1)))
Y2 <- Al * U2 + V2 + as.vector(rnorm(n2, mean = 0, sd = sqrt(tau2)))


########################################################################
##########  Recover parameters with full conditionals defined
########################################################################

#######################Codes below are old codes that need update
# Simulation for U1: 
set.seed(123)
U2 <- as.vector(t(chol(sig1 * S22)) %*% rnorm(n2))
S1 <- S11 - S12 %*% solve(S22) %*% S21
A1 <- S12 %*% solve(S22)
S1inv <- solve(S1)

SigmaU1 <- solve(1/tau1 * diag(1, n1) + 1/sig1 * S1inv)
meanU1 <- SigmaU1 %*% (1/tau1 * Y1 + 1/sig1 * S1inv %*% A1 %*% U2)

U1_sim <- as.vector(t(chol(SigmaU1)) %*% rnorm(n1)) + meanU1


# Simulation for U2:
# A_l = rho in this setting
set.seed(123)
U1 <- as.vector(t(chol(sig1 * S11)) %*% rnorm(n1))
V2 <- as.vector(t(chol(sig2 * S22)) %*% rnorm(n2))
S2 <- S22 - S21 %*% solve(S11) %*% S12
A2 <- S21 %*% solve(S11)
S2inv <- solve(S2)

SigmaU2 <- solve(rho^2/tau2 * diag(1, n2) + 1/sig1 * S2inv)
meanU2 <- SigmaU2 %*% (1/tau2 * rho * Y2 - rho * V2 + 1/sig1 * S2inv %*% A2 %*% U1)

U2_sim <- as.vector(t(chol(SigmaU2)) %*% rnorm(n2)) + meanU2
hist(U2_sim)

# Simulation for V2:
set.seed(123)
SigmaV2 <- solve(1/sig2 * solve(S22) + 1/tau2 * diag(1, n2))
meanV2 <- SigmaV2 %*% (1/tau2 * Y2)

V2_sim <- as.vector(t(chol(SigmaV2)) %*% rnorm(n2)) + meanV2
hist(V2_sim)





















