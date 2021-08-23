rm(list=ls())

#######################################################
#
# INPUTS
#
# Y1         := n1 x nt-vector of observations
# Y2         := n2 x nt-vector of observations
# s1         := n1x2 matrix of spatial coordinates
# s2         := n2x2 matrix of spatial coordinates
# mean_range := log(spatial range) ~ N(mean_range,sd_range)
# sd_range
# mean_var   := log(variance parameters) ~ N(mean_var,sd_var)
# sd_var
# mean_rho   := rho ~ N(mean_rho,sd_rho)
# sd_rho
#
# iters      := number of MCMC iterations
# burn       := length of burn-in
#
# Model
#
# Y1[i,t]    = beta1 + Z1[i,t] + e1[i,t]
# Y2[i,t]    = beta2 + rho*Z1[i,t] + Z2[i,t] + e2[i,t]
# Var(Z) = sigU
# Var(Z) = sigV
# Var(e1) = tau1
# Var(e2) = tau2
# Cor(Z1(s1),Z1(s2)) = exp(-d12/range1)
# Cor(Z2(s1),Z2(s2)) = exp(-d12/range2)
#
#######################################################


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


# The main MCMC function
LMC <- function(Y1,Y2,s1,s2,
                mean_range=0,sd_range=1,
                mean_var=0,sd_var=1,
                mean_rho=0,sd_rho=10,
                iters=5000,burn=1000,thin=1,update=10){

   library(fields) 
   library(geoR)

   tick <- proc.time()[3]

   # Bookkeeping
    n1       <- nrow(Y1)
    n2       <- nrow(Y2)
    nt       <- ncol(Y1)
    m1       <- is.na(Y1)
    m2       <- is.na(Y2)
    d        <- as.matrix(dist(rbind(s1,s2)))
    const    <- 1/nt

   # Initial values

    beta1 <- mean(Y1,na.rm=TRUE)
    beta2 <- mean(Y2,na.rm=TRUE)
    Y1[m1] <- beta1
    Y2[m2] <- beta2
    U1     <- matrix(0,n1,nt)
    U2     <- matrix(0,n2,nt)
    V2     <- matrix(0,n2,nt)
    rhoU   <- exp(mean_rho)
    rhoV   <- exp(mean_rho)
    taue1  <- var(as.vector(Y1))
    taue2  <- var(as.vector(Y2))
    tauU   <- rep(taue1,nt)
    tauV   <- rep(taue1,nt)
    A      <- rep(1,nt)
    Z1     <- matrix(0,n1,nt)
    Z2     <- matrix(0,n1,nt)
    Ys1    <- matrix(0,n1,nt)
    Ys2    <- matrix(0,n1,nt)
    eigU   <- invU(d,n1,n2,rhoU)
    eigV   <- invV(d,rhoV)


   # Keep track of stuff

    keep_theta  <- matrix(0,iters,6)
    colnames(keep_theta) <- c("rhoU","rhoU","tauU","tauV","taue1","taue2")

    att <- acc <- MH <- rep(0.5,q)

   # GO!!!

   for(iter in 1:iters){


     for(ttt in 1:thin){

      ##############################################:
      ####     Transform to spatial land       #####:
      ##############################################:

        for(i in 1:n1){Z1[i,] <- fft_real(U1[i,],inverse=TRUE)}
        for(i in 1:n2){Z2[i,] <- fft_real(A*U2[i,]+V2[i,],inverse=TRUE)}

      ##############################################:
      ####  IMPUTE MISSING DATA (real space)   #####:
      ##############################################:

        Y1[m1] <- rnorm(sum(m1),beta1+Z1[m1],1/sqrt(taue1))
        Y2[m2] <- rnorm(sum(m2),beta2+Z2[m2],1/sqrt(taue2))
 
      ##############################################:
      ####      MEAN PARAMETERS (real space)   #####:
      ##############################################:

        VVV   <- taue1*n1*nt + 0.01
        MMM   <- taue1*sum(Y1-Z1) 
        beta1 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
        
        VVV   <- taue2*n2*nt + 0.01
        MMM   <- taue2*sum(Y2-Z2) 
        beta2 <- rnorm(1,MMM/VVV,1/sqrt(VVV))

      ##############################################:
      ####     ERROR VARIANCES (real space)    #####:
      ##############################################:

        taue1 <- rgamma(1,n1*nt/2+0.01,sum((Y1-beta1-Z1)^2)/2+0.01)
        taue2 <- rgamma(1,n2*nt/2+0.01,sum((Y2-beta2-Z2)^2)/2+0.01)

      ##############################################:
      ####     Transform to spectral land      #####:
      ##############################################:

        for(i in 1:n1){Ys1[i,] <- fft_real(Y1[i]-beta1)}
        for(i in 1:n2){Ys2[i,] <- fft_real(Y2[i]-beta2)}
        taus1 <- const*taue1
        taus2 <- const*taue2

      ##############################################:
      ####      LMC TERMS (spectral space)     #####:
      ##############################################:

# Here is where we would put the complicated full conditions that are
# derived in overleaf

      ###################################################:
      ####  COVARIANCE PARAMETERS (spectral domain) #####:
      ###################################################:

       # tauU ... bV are fixed here but should be updated
       tauA <- muA  <- 0
       aU   <-  bU <- aV <- bV <- 0.1
       for(t in 1:nt){
         VVV     <- taus2*sum(U2[,t]^2)+tauA
         MMM     <- taus2*sum((Ys2[,t]-V2[,t])*U2[,t])+muA*tauA
         A[t]    <- rnorm(1,MMM/VVV,1/sqrt(VVV))   

         U       <- c(U1[,t],U2[,t])
         tauU[t] <- rgamma(1,(n1+n2)/2+aU,t(U)%*%E$inv%*%U/2+bU)
         tauV[t] <- rgamma(1,n2/2+aV,t(V2[,t])%*%E$inv2%*%V2[,t]/2+bV)
       }


     } # end thin

      ##############################################:
      #####        KEEP TRACK OF STUFF       #######:
      ##############################################:

#       keep_theta[iter,] <- ...

      ##############################################:
      #####       PLOT RESULTS SO FAR        #######:
      ##############################################:

      if(iter%%update==0){
        par(mfrow=c(3,3))
        for(j in 1:2){
         plot(keep_beta[1:iter,j],type="l",
              xlab="MCMC iteration",ylab="Sample",
              main=colnames(keep_beta)[j])
        }
        for(j in 1:q){
         plot(keep_theta[1:iter,j],type="l",
              xlab="MCMC iteration",ylab="Sample",
              main=colnames(keep_theta)[j])
        }
      }


   }   
   tock <- proc.time()[3]
 out <- list(beta=keep_beta,theta=keep_theta,
             acc_rate=acc/att,time=tock-tick)
return(out)}


# Generate a fake dataset

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
 beta   <- c(4,-2)
 Int1   <- ifelse(type==1,1,0)
 Int2   <- ifelse(type==2,1,0)
 X      <- cbind(Int1,Int2)
 Sig    <- cov.LMC(as.matrix(dist(s)),type,theta)
 Y      <- as.vector(X%*%beta + t(chol(Sig$S))%*%rnorm(n))

# Fit the model
 
 fit    <- LMC(Y,X,s,type)

# Plot the results

 par(mfrow=c(3,3))
 for(j in 1:2){
   plot(fit$beta[,j],type="l",
        xlab="MCMC iteration",ylab="Sample",
        main=colnames(fit$beta)[j])
   abline(beta[j],0,col=2,lwd=2) 
 }
 for(j in 1:length(theta)){
   plot(fit$theta[,j],type="l",
        xlab="MCMC iteration",ylab="Sample",
        main=colnames(fit$theta)[j])
   abline(theta[j],0,col=2,lwd=2) 
 }

