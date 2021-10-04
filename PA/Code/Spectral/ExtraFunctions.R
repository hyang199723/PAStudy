##############################################################################
                      ##### Extra Functions ######
###############################################################################
                              
# Compute exponential correlation function
# d: distance matrix between coords
# range: range parameter
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
  S[t==2,t==2] <- sig1*exp(-d[t==2,t==2]/range1) # Should it be range1 or range2?
  return(S)
}
covUr <- function(d,t,range) {
  S <- diag(0, length(t))
  S[t==1,t==1] <- exp(-d[t==1,t==1]/range)
  S[t==1,t==2] <- exp(-d[t==1,t==2]/range)
  S[t==2,t==1] <- exp(-d[t==2,t==1]/range)
  S[t==2,t==2] <- exp(-d[t==2,t==2]/range) 
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


covV2 <- function(d,t,lrange2) {
  # We only need type = 2
  range2 <- exp(lrange2)
  d <- d[t == 2, t== 2]
  t <- t[t == 2]
  S <- diag(0, length(t))
  S <- exp(-d/range2)
  return(S)
}

# Function to compute Sigma11, Sigma12, Sigma21, and Sigma22
# d: distance matrix between coords
# t: type, either 1 or 2
# theta: parameters
# # theta  <- c(rho,log(range1),log(range2),log(tau1),log(tau2),log(sig1),log(sig2))
covIndividual <- function(d, t, theta2) {
  range1 <- exp(theta2)
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
  return(out)}


# The log likelihood function
log_like <- function(Y,SIG){
  -0.5*SIG$logdet - 0.5*sum(Y*(SIG$Sinv%*%Y))
}

invU <- function(d,n1,n2,rho){
  S   <- exp(-d/rho)
  E   <- eigen(S)
  G   <- E$vectors
  D   <- E$values
  S1  <- S[1:n1,1:n1]
  S2  <- S[1:n2+n1,1:n2+n1]
  S12 <- S[1:n1,1:n2+n1]
  Q   <- G%*%diag(1/D)%*%t(G)
  Q1  <- Q[1:n1,1:n1]
  Q2  <- Q[1:n2+n1,1:n2+n1]
  A12 <- S12%*%solve(S2)
  A21 <- t(S12)%*%solve(S1)
  out <- list(G=G,D=D,Q=Q,Q1=Q1,Q2=Q2,A12=A12,A21=A21)
  return(out)}


invV <- function(d,rho){
  S   <- exp(-d/rho)
  E   <- eigen(S)
  G   <- E$vectors
  D   <- E$values
  Q   <- G%*%diag(1/D)%*%t(G)
  out <- list(G=G,D=D,Q=Q)
  return(out)}
