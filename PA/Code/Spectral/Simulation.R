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
  range2 <- exp(theta[3])
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
  range2 <- exp(theta[3])
  S <- diag(0, length(t))
  S[t==1,t==1] <- exp(-d[t==1,t==1]/range1)
  S[t==1,t==2] <- exp(-d[t==1,t==2]/range1)
  S[t==2,t==1] <- exp(-d[t==2,t==1]/range1)
  S[t==2,t==2] <- exp(-d[t==2,t==2]/range1)
  S11 <- S[t==1,t==1]
  S12 <- S[t==1,t==2]
  S21 <- S[t==2,t==1]
  S22 <- S[t==2,t==2]
  return(list(S11=S11, S12=S12, S21=S21, S22=S22))
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
tau1   <- 0.1
tau2   <- 0.3
set.seed(125); sig1 <- 1#rinvgamma(1, 2, 1) # ~invgamma(2,1)
set.seed(124); sig2 <- 2#rinvgamma(1, 2, 2) # ~invgamma(2,2)
theta  <- c(rho,log(range1),log(range2),
            log(tau1),log(tau2),log(sig1),log(sig2))
n1 <- sum(type == 1)
n2 <- sum(type == 2)

# Covariance matrix for U1, U2
Su <- covU(dist,type,theta)
# Covariance matrix for V2
Sv <- covV(dist,type,theta)

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
set.seed(12)
Al <- rtruncnorm(1, a=0, b=+Inf, mean=2, sd=2)

# Generate Y1 and Y2
set.seed(123)
Y1 <- U1 + as.vector(rnorm(n1, mean = 0, sd = sqrt(tau1)))
Y2 <- Al * U2 + V2 + as.vector(rnorm(n2, mean = 0, sd = sqrt(tau2)))


########################################################################
##########  Recover parameters with full conditionals defined
########################################################################
iter = 3000
S <- covIndividual(dist, type, theta)
S11 <- S$S11; S12 <- S$S12; S21 <- S$S21; S22 <- S$S22
# Simulation for U1: 
S1 <- S11 - S12 %*% solve(S22) %*% S21
S1inv <- solve(S1)
A1 <- S12 %*% solve(S22)
sigmaU1 <- solve(1/tau1 * diag(1, n1) + 1/sig1 * S1inv)
meanU1 <- sigmaU1 %*% (1/tau1 * Y1 + 1/sig1 * S1inv %*% A1 %*% U2)

U1_sim_all <- matrix(data = 0, nrow = n1, ncol = iter)
for (i in 1:iter) {
  U1_sim <- as.vector(t(chol(sigmaU1)) %*% rnorm(n1)) + meanU1
  U1_sim_all[,i] = U1_sim
}
U1_sim <- as.vector(apply(U1_sim_all, 1, mean))
diff <- (t(U1_sim - U1) %*% (U1_sim - U1)) / (t(U1) %*% U1) # 4% difference
angel <- t(U1_sim) %*% U1 / sqrt((t(U1) %*% U1)) / sqrt((t(U1_sim) %*% U1_sim)) # Cos() = 0.9816

# Simulation for U2:
S2 <- S22 - S21 %*% solve(S11) %*% S12
A2 <- S21 %*% solve(S11)
S2inv <- solve(S2)

sigmaU2 <- solve(Al^2/tau2 * diag(1, n2) + 1/sig1 * S2inv)
meanU2 <- sigmaU2 %*% (1/tau2 * Al * Y2 - 1/tau2 * Al * V2 + 1/sig1 * S2inv %*% A2 %*% U1)

U2_sim_all <- matrix(data = 0, nrow = n2, ncol = iter)
for (i in 1:iter) {
  U2_sim <- as.vector(t(chol(sigmaU2)) %*% rnorm(n2)) + meanU2
  U2_sim_all[,i] = U2_sim
}
U2_sim <- as.vector(apply(U2_sim_all, 1, mean))
diff <- (t(U2_sim - U2) %*% (U2_sim - U2)) / (t(U2) %*% U2) # 58% difference
angel <- t(U2_sim) %*% U2 / sqrt((t(U2) %*% U2)) / sqrt((t(U2_sim) %*% U2_sim)) # Cos() = 0.5

# Simulation for V2:
#Sv = sig2 * Sigma_{v22}
sigmaV2 <- solve(solve(Sv) + 1/tau2 * diag(1, n2))
meanV2 <- sigmaV2 %*% (1/tau2 * Y2 - 1/tau2 * Al * U2)

V2_sim_all <- matrix(data = 0, nrow = n2, ncol = iter)
for (i in 1:iter) {
  V2_sim <- as.vector(t(chol(sigmaV2)) %*% rnorm(n2)) + meanV2
  V2_sim_all[,i] = V2_sim
}
V2_sim <- as.vector(apply(V2_sim_all, 1, mean))
diff <- (t(V2_sim - V2) %*% (V2_sim - V2)) / (t(V2) %*% V2) # 0.349 difference
angel <- t(V2_sim) %*% V2 / sqrt((t(V2) %*% V2)) / sqrt((t(V2_sim) %*% V2_sim)) # Cos() = 0.95



# Simulation for Al; Al ~ TN(2,2); Al = 0.87
sigmaAl <- solve(1/tau2 * t(U2) %*% U2 + 1/4)
meanAl <- sigmaAl %*% (1/tau2 + t(Y2) %*% U2 - 1/tau2 * t(V2) %*% U2 + 2/4) # Mean is 0.879
Al_sim <- rtruncnorm(3000, a=0, b=+Inf, mean=meanAl, sd=sigmaAl)
hist(Al_sim)
mean(Al_sim) # 1.48

# Simulation for sig1; sig1~IG(2,1)
a <- (n2+n1)/2 + 2
b <- (t(U) %*% solve(Su/sig1) %*% U) / 2 + 1
sig2_sim <- rinvgamma(3000, a, b)
a/(b-1) # 1.149


# Simulation for sig2; Sig2~IG(2,2)
a <- n2/2 + 2
b <- (t(V2) %*% solve(Sv/sig2) %*% V2) / 2 + 2
# Mean = 3.14






