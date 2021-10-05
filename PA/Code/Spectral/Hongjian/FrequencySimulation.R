rm(list = ls())
library(mvtnorm)
library(MASS)
library(truncnorm)
library(invgamma)
setwd("/Users/hongjianyang/Research/PAStudy/PA/")
raw_frm <- read.csv("Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv", header = TRUE)
raw_pa <- read.csv("Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv", header = TRUE)
raw_frm$Timestamp <- as.POSIXct(raw_frm$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
raw_pa$Timestamp <- as.POSIXct(raw_pa$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
# Select three timestamps
start = as.POSIXct('2020-01-30 12:00:00')
end = as.POSIXct('2020-01-30 14:00:00')
frm <- subset(raw_frm, (Timestamp >= start) & (Timestamp <= end))
pa <- subset(raw_pa, (Timestamp >= start) & (Timestamp <= end))
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


fft_real <- function(dat,inverse=FALSE){
  if(!inverse){
    x  <- dat
    n  <- length(x)
    n2 <- floor(n/2)
    y  <- fft(x,inverse=FALSE)
    if(n%%2==0){
      X1     <- Re(y)[1:(n2+1)]
      X2     <- Im(y)[2:(n2)]
    }
    if(n%%2!=0){
      X1     <- Re(y)[1:(n2+1)]
      X2     <- Im(y)[2:(n2+1)]
    }
    out <- c(X1,X2)
  }
  if(inverse){
    X  <- dat
    n  <- length(X)
    n2 <- floor(n/2)
    if(n%%2==0){
      Y1    <- c(X[1:(n2+1)],X[n2:2])
      Y2    <- c(0,X[(n2+2):n],0,-X[n:(n2+2)])
    }
    if(n%%2!=0){
      Y1    <- c(X[1:(n2+1)],X[(n2+1):2])
      Y2    <- c(0,X[(n2+2):n],-X[n:(n2+2)])
    }
    y   <- complex(n, real = Y1, imaginary = Y2)
    out <- Re(fft(y/n,inverse=TRUE))
  }
  return(out)
}

########################################################################
##########      Simulate data on 5 timestamps
########################################################################
# Total number of observations
n      <- 250
nt <- 5
# data matrix
dat <- matrix(0, nrow = n, ncol = nt)
# 70% PA(Type 2), 30% AQS (Type 1)
type   <- rbinom(n,1,0.7)+1
# Ordered, easier for simulation
type <- type[order(type)]
# Generate coords
s      <- cbind(runif(n),runif(n))
# Distance matrix
dist <- as.matrix(dist(s))
# Parameter vectors
# \sigma_{ul} = \sig1, #\sigma_{vl} = \sig2
rho_vec <- rnorm(nt)
range1_vec <- runif(nt)
range2_vec <- runif(nt, 0, 2)
tau1_vec <- runif(nt)
tau2_vec <- runif(nt, 0, 2)
sig1_vec <- runif(nt)
sig2_vec <- runif(nt, 0, 2)
Al_vec <- runif(nt, 0, 2)
n1 <- sum(type == 1)
n2 <- sum(type == 2)
for (i in 1:nt) {
  rho    <- rho_vec[i]
  range1 <- range1_vec[i]
  range2 <- range2_vec[i]
  tau1   <- tau1_vec[i]
  tau2   <- tau2_vec[i]
  sig1 <- sig1_vec[i]
  sig2 <- sig2_vec[i]
  # parameter vector
  theta  <- c(rho,log(range1),log(range2),
              log(tau1),log(tau2),log(sig1),log(sig2))
  
  # Covariance matrix for U1, U2
  Su <- covU(dist,type,theta)
  Su1 <- covIndividual(dist,type,theta)
  # Covariance matrix for V2
  Sv <- covV(dist,type,theta)
  Sv1 <- Sv/sig2
  # Generate U1, U2, V2, Al
  # Generate U1, U2
  # The first n1 values are U1, remaining are U2
  U <- as.vector(t(chol(Su1$S)) %*% rnorm(n1 + n2, mean = 0, sd = sqrt(sig1)))
  U1 <- U[1:n1]
  U2 <- U[n1+1:n2]
  # Generate V2
  V2 <- as.vector(t(chol(Sv1)) %*% rnorm(n2, mean = 0, sd = sqrt(sig2)))
  # Generate A_{l}
  Al <- Al_vec[i]
  # Generate Y1 and Y2
  Y1 <- as.vector(U1 + as.vector(rnorm(n1, mean = 0, sd = sqrt(tau1))))
  Y2 <- as.vector(Al * U2 + V2 + as.vector(rnorm(n2, mean = 0, sd = sqrt(tau2))))
  dat[1:n1, i] = Y1
  dat[(n1+1):(n1+n2), i] = Y2
}

########################################################################
##########      Convert to frequency domain
########################################################################
# Bookkeeping
iter <- 3000
Ys1 <- matrix(0, nrow = n1, ncol = nt)
Ys2 <- matrix(0, nrow = n2, ncol = nt)
for (i in 1:nt) {
  Ys1[, i] = fft_real[dat[1:n1, i]]
  Ys2[, i] = fft_real[dat[(n1+1):(n1+n2), i]]
}



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
  
  # Update range1 and range2
  iters=5000
  nchains=2
  
  keep.range1=matrix(NA,iters,nchains)
  keep.range2=matrix(NA,iters,nchains)
  
  # Initial values
  
  Ru = as.vector(u)/sqrt(sigmau)
  Rv = as.vector(v2)/sqrt(sigmav)
  
  for (j in 1:nchains)
  {
    
    lrange1 = 0
    lrange2 = 0
    for(i in 1:iters){
      # range1
      Ms=exp_corr(du12,range=exp(lrange1))
      curll = dmvnorm(Ru,rep(0,sum(n)),Ms,log=TRUE)
      canrange1 = rnorm(1,lrange1,0.5)
      canM = exp_corr(du12,range=exp(canrange1))
      canll = dmvnorm(Ru,rep(0,sum(n)),canM,log=TRUE)
      
      MH1 <- canll-curll+dnorm(canrange1,log=TRUE)-dnorm(lrange1,log=TRUE)
      
      if (log(runif(1))<MH1)
      {
        lrange1=canrange1
      }
      keep.range1[i,j]  <-exp(lrange1)
      
      # range2
      Ss=exp_corr(dv2,range = exp(lrange2))
      curll2 = dmvnorm(Rv,rep(0,n[2]),Ss,log=TRUE)
      canrange2 = rnorm(1,lrange2,0.5)
      canS = exp_corr(dv2,range=exp(canrange2))
      canll2 = dmvnorm(Rv,rep(0,n[2]),canS,log=TRUE)
      
      MH2 <- canll2-curll2+dnorm(canrange2,log=TRUE)-dnorm(lrange2,log=TRUE)
      
      if (log(runif(1))<MH2)
      {
        lrange2=canrange2
      }
      keep.range2[i,j]  <-exp(lrange2)
      
      
    }
}
