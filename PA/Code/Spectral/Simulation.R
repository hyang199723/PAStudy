# Function to compute the MV spatial covariance matrix
cov.LMC <- function(d,t,theta){
  # t[i] = 1 means AQS; t[i] = 2 mean PA
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


set.seed(123)
n      <- 250
type   <- rbinom(n,1,0.3)+1
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
