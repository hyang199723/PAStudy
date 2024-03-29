# All LMC functions
# Last updated: 09/23/2022
library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)
library(mvtnorm)
library(truncnorm)


exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
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
  S2  <- S[(n1+1):(n2+n1),(n1+1):(n2+n1)]
  S12 <- S[1:n1,(n1+1):(n2+n1)]
  S21 <- S[(n1+1):(n2+n1),1:n1]
  Q   <- G%*%diag(1/D)%*%t(G)
  Q1  <- Q[1:n1,1:n1]
  Q2  <- Q[1:n2+n1,1:n2+n1]
  A12 <- S12%*%solve(S2)
  A21 <- t(S12)%*%solve(S1)
  out <- list(G=G,D=D,Q=Q,Q1=Q1,Q2=Q2,A12=A12,A21=A21,S11=S1,S12=S12,S22=S2,S21=S21)
  return(out)
}


invV <- function(d,rho){
  S   <- exp(-d/rho)
  E   <- eigen(S)
  G   <- E$vectors
  D   <- E$values
  Q   <- G%*%diag(1/D)%*%t(G)
  out <- list(G=G,D=D,Q=Q)
  return(out)
}

# Input specification:
# Y1: Type 1 response value; Matrix with (number of observation) * (number of time steps)
#     Y1 can have missing values
# Y2: Type 2 response value; Matrix with (number of observation) * (number of time steps)
#     Y2 can have missing values
# s1: locations of Type 1 response
# s2: locations of Type 2 response
# sp1: prediction locations for Type 1 response
LMC_fit=function(Y1,Y2, s1,s2,sp1=NULL,
                 mean_range=0, sd_range=1, mean_var=0, sd_var=1, mean_rho=0,
                 sd_rho=10, iters=3000, burn=1000, thin=1, update=10)
{
  n1       <- nrow(Y1)
  n2       <- nrow(Y2)
  nt       <- ncol(Y1)
  m1       <- is.na(Y1)
  m2       <- is.na(Y2)
  d        <- as.matrix(dist(rbind(s1,s2)))
  dv2 = as.matrix(dist(s2))
  const    <- nt/2
  # constant for sigmaU-V priors 
  a1 <- (n2+n1)/2 + 1
  a2 <- n2/2 + 1  
  
  # Keep track of stuff
  
  # All parameters: "rangeu","rangev","sigmaU","sigmaV","taue1","taue2","A",'beta1','beta2'
  
  # Parameters that are different for each spectrum: sigmaU, sigmaV, A, beta1, beta2
  # Dimension of data structure: nt * iterations
  keep.sigmaU = array(0, dim = c(nt, iters,thin))
  keep.sigmaV = array(0, dim = c(nt, iters,thin))
  keep.A = array(0, dim = c(nt, iters,thin))
  
  # Parameters that stay the same: rangeU, rangeV, taus1, taus2
  # Dimension of data structure: 1 * iterations
  keep.rangeU = matrix(0, iters,thin)
  keep.rangeV = matrix(0, iters,thin)
  # keep.taus1 = matrix(0, iters,thin)
  # keep.taus2 = matrix(0, iters,thin)
  keep.taue1 = matrix(0, iters,thin)
  keep.taue2 = matrix(0, iters,thin)
  # Vector parameters; Dimension of data structure: #stations * nt * iters
  # keep.u1= array(0,dim=c(n1,nt,iters,thin))
  # keep.u2= array(0,dim=c(n2,nt,iters,thin))
  # keep.v2= array(0,dim=c(n2,nt,iters,thin))
  #keep.Y1.M= array(0,dim=c(n1,nt,iters,thin))
  #keep.Y2.M= array(0,dim=c(n2,nt,iters,thin))
  
  # Save Mean and Variance of U and V 
  E.u1=matrix(0,ncol=nt,nrow=n1)
  E.u2=matrix(0,ncol=nt,nrow=n2)
  E.v=matrix(0,ncol=nt,nrow=n2)

  E2.u1=matrix(0,ncol=nt,nrow=n1)
  E2.u2=matrix(0,ncol=nt,nrow=n2)
  E2.v=matrix(0,ncol=nt,nrow=n2)
  
  E.corrY1Y2=numeric(nt)
  E2.corrY1Y2=numeric(nt)
  
  ## get info for the ranges priors 
  priorR_mn1 <- log(max(d)) - 1.5
  priorR_sd1 <- .5
  
  priorR_mn2 <- log(max(dv2)) - 1.5
  priorR_sd2 <- .5
  
  predictions <- !is.null(sp1)
  if(predictions){
    np1 = nrow(sp1)
    all.d=as.matrix(dist(rbind(s1,s2,sp1))) 
    keep.Y1.P= array(0,dim=c(np1,nt,iters,thin))
  }
  
  aru=rep(0,thin)
  arv=rep(0,thin)
  
  # start MCMC
  
  for(ttt in 1:thin){
    # Mean imputation
    beta1 <- mean(colMeans(Y1,na.rm=TRUE),na.rm=TRUE)
    beta2 <- mean(colMeans(Y2,na.rm=TRUE),na.rm=TRUE)
    Y1[m1] <- beta1
    Y2[m2] <- beta2
    # create initial vectors 
    U1     <- matrix(0,n1,nt)
    U2     <- matrix(0,n2,nt)
    V2     <- matrix(0,n2,nt)
    Z1     <- matrix(0,n1,nt)
    Z2     <- matrix(0,n2,nt)
    Ys1    <- matrix(0,n1,nt)
    Ys2    <- matrix(0,n2,nt)
    
    # Initial parameter values
    rangeU  <- 0.4
    rangeV   <- 0.2
    lrangeU  <- log(rangeU)
    lrangeV   <-log(rangeV)
    
    taue1  <- var(as.vector(Y1))
    taue2  <- var(as.vector(Y2))
    sigmaU   <- rep(taue1,nt)
    sigmaV   <- rep(taue2,nt)
    A      <- rep(0,nt)
    
    if(predictions){
      # set initial values
      U1p = matrix(0,np1,nt)
      Z1p = matrix(0,np1,nt)
      Ys1.pred = matrix(0,np1,nt)
      Y1.pred = matrix(beta1,np1,nt)
    }
    
    for(iter in 1:iters){
      
      
      ##############################################:
      ####     Transform to spatial land       #####:
      ##############################################:
      
      for(i in 1:n1){Z1[i,] <- fft_real(U1[i,],inverse=TRUE)}
      for(i in 1:n2){Z2[i,] <- fft_real(A*U2[i,]+V2[i,],inverse=TRUE)}
      
      ##############################################:
      ####  IMPUTE MISSING DATA (real space)   #####:
      ##############################################:
      
      Y1[m1] <- rnorm(sum(m1),beta1+Z1[m1],sqrt(taue1)) 
      Y2[m2] <- rnorm(sum(m2),beta2+Z2[m2],sqrt(taue2)) 
      
      ##############################################:
      ####      MEAN PARAMETERS (real space)   #####:
      ##############################################:
      
      # full conditional for beta1 and beta2 ...
      #VVV   <- (taue1*n1*nt + 0.01)
      #MMM   <- taue1*sum(Y1-Z1) 
      #beta1 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      #beta1=beta.1
      
      # VVV   <- (taue2*n2*nt + 0.01)
      # MMM   <- taue2*sum(Y2-Z2) 
      #beta2 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      #beta2=beta.2
      
      ##############################################:
      ####     ERROR VARIANCES (real space)    #####:
      ##############################################:
      # full conditionals for taue1 and taue2 
      taue1 <- 1/rgamma(1,n1*nt/2+1,sum((Y1-beta1-Z1)^2)/2+1)
      taue2 <- 1/rgamma(1,n2*nt/2+1,sum((Y2-beta2-Z2)^2)/2+1)
      
      ##############################################:
      ####     Transform to spectral land      #####:
      ##############################################:
      
      for(i in 1:n1){Ys1[i,] <- fft_real(as.numeric(Y1[i,] - beta1))}
      for(i in 1:n2){Ys2[i,] <- fft_real(as.numeric(Y2[i,] - beta2))}
      taus1 <- nt / 2 *taue1
      taus2 <- nt / 2 *taue2
      
      ##############################################:
      ####      LMC TERMS (spectral space)     #####:
      ##############################################:
      
      eigU = invU(d,n1,n2,rangeU)
      S11 = eigU$S11; S12 <- eigU$S12; S21 <- eigU$S21; S22 <- eigU$S22
      S1 = S11 - eigU$A12 %*% S21
      ES1=eigen(S1)
      S1_G   = ES1$vectors
      S1_D   = ES1$values
      S1inv = S1_G%*%diag(1/S1_D)%*%t(S1_G)
      A1 = S12 %*% solve(S22)
      S1invA1 = S1inv %*% A1
      # S1inv = G * diag(1/D) * G'
      # inverse = G * [1/taus1 + 1/sigmaU * diag(1/D)]^{-1} * G'
      
      S11inv=solve(S11)
      S2 = S22 - S21 %*% S11inv %*% S12
      A2 = S21 %*% S11inv
      ES2 = eigen(S2)
      S2_G = ES2$vectors
      S2_D = ES2$values
      S2inv = S2_G%*%diag(1/S2_D)%*%t(S2_G)
      S2invA2 = S2inv %*% A2
      
      # G: eigenvectors
      # D: eigenvalues
      # Q: inverse matrix
      eigV  = invV(dv2,rangeV)
      V_G = eigV$G
      V_D = eigV$D
      
      S2Vinv= invV(dv2,rangeV)$Q
      
      GYs1 = t(S1_G)%*%Ys1
      GSU2 = t(S1_G)%*%S1invA1%*%U2
      GVYs2 = t(V_G)%*%Ys2
      GYs2 = t(S2_G)%*%Ys2
      GSU1 = t(S2_G)%*%S2invA2%*%U1
      GSV2 = t(S2_G)%*%V2
      GVU2 <- t(V_G)%*%U2
      
      # Gibbs sampling
      for (r in 1:nt) # for each spectral
      {
        # # Sample U1
        newDiag = 1/taus1 + 1/sigmaU[r] * 1/S1_D
        sigmaU1 <- S1_G %*% diag(1/newDiag) %*% t(S1_G)
        meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1invA1 %*% U2[,r])
        U1[,r] = meanU1+S1_G%*%rnorm(n1, 0, 1/sqrt(newDiag))
        # meanU1 = GYs1[,r]/taus1 + GSU2[,r]/sigmaU[r]
        # U1[,r] = S1_G%*%(meanU1/newDiag + rnorm(n1, 0, 1/sqrt(newDiag)))
        
        # # Sample U2
        newDiag = A[r]^2/taus2 + 1/sigmaU[r] * 1/S2_D
        sigmaU2 <- S2_G %*% diag(1/newDiag) %*% t(S2_G)
        meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 *A[r] * V2[,r] +
                                 1/sigmaU[r] * S2invA2 %*% U1[,r])
        U2[,r] = meanU2+S2_G%*%rnorm(n2,0,1/sqrt(newDiag))
        # meanU2 <- (GYs2[,r]*A[r]*1/taus2-GSV2[,r]*A[r]/taus2+GSU1[,r]*1/sigmaU[r])
        # U2[,r] <- S2_G%*%(meanU2/newDiag+rnorm(n2, 0, 1/sqrt(newDiag)))
        
        # # Sample V2
        newDiag = 1/sigmaV[r] * 1/V_D + 1/taus2
        sigmaV2 = V_G %*% diag(1/newDiag) %*% t(V_G)
        meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
        V2[,r] =meanV2+V_G%*%rnorm(n2,0,1/sqrt(newDiag))
        # meanV2 <- (1/taus2 * GVYs2[,r] - 1/taus2 * A[r] * GVU2[,r])
        # V2[,r] <- V_G%*%(meanV2/newDiag+rnorm(n2, 0, 1/sqrt(newDiag)))
        
        # # Sample Al
        sigmaAl <- solve(1/taus2 * sum(U2[,r]^2) + 1/5)
        meanAl <- sigmaAl %*% (1/taus2 *sum(Ys2[,r]*U2[,r]) - 1/taus2 * sum(V2[,r]*U2[,r]) + 0.8/5)
        A[r] <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl, sd=sqrt(sigmaAl))
        
        # # Sample sig1
        U_sim <- as.vector(append(U1[,r], U2[,r]))
        b1 <- (t(U_sim) %*% eigU$Q %*% U_sim) / 2 + 1
        sigmaU[r] = 1/rgamma(1, a1, b1)
        
        # # # Sample sig2
        b2 <- (t(V2[,r]) %*% eigV$Q %*% V2[,r]) / 2 + 1
        sigmaV[r] = 1/rgamma(1, a2, b2)
        
      }
      
      ###############################################
      ##       Metropolis H: Range parameters      ##
      ###############################################
      # Sweep operation to get rid of variance
      Ru=sweep(rbind(U1,U2),2,FUN='/',sqrt(sigmaU))
      Rv=sweep(V2,2,FUN='/',sqrt(sigmaV))
      
      # range1
      Ms=exp_corr(d,range=exp(lrangeU))
      #curll = sum(apply(Ru,2,dmvnorm,mean=rep(0,n1+n2),sigma=Ms,log=TRUE))
      curll = sum(dmvnorm(t(Ru), mean=rep(0,n1+n2),sigma=Ms,log=TRUE))
      prior_curll=dnorm(lrangeU,mean=priorR_mn1,sd=priorR_sd1,log=TRUE)
      
      canrange1 = rnorm(1,lrangeU,0.1)
      canM = exp_corr(d,range=exp(canrange1))
      #canll = sum(apply(Ru,2,dmvnorm,mean=rep(0,n1+n2),sigma=canM,log=TRUE))
      canll = sum(dmvnorm(t(Ru), mean=rep(0,n1+n2),sigma=canM,log=TRUE))
      prior_canll=dnorm(canrange1,mean=priorR_mn1,sd=priorR_sd1,log=TRUE)
      
      MH1 <- canll-curll+prior_canll-prior_curll#+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
      
      if (log(runif(1))<MH1)
      {
        lrangeU=canrange1
        aru[ttt]=aru[ttt]+1
      }
      
      # range2
      Ss=exp_corr(dv2,range = exp(lrangeV))
      #curll2 = sum(apply(Rv,2,dmvnorm,mean=rep(0,n2),sigma=Ss,log=TRUE))
      curll2 = sum(dmvnorm(t(Rv), mean = rep(0,n2),sigma=Ss,log=TRUE ))
      prior_curll2=dnorm(lrangeV,mean=priorR_mn2,sd=priorR_sd2,log=TRUE)
      
      canrange2 = rnorm(1,lrangeV,0.1)
      canS = exp_corr(dv2,range=exp(canrange2))
      #canll2 = dmvnorm(Rv,rep(0,n2),canS,log=TRUE)
      canll2 = sum(dmvnorm(t(Rv), mean = rep(0,n2),sigma=canS,log=TRUE ))
      prior_canll2=dnorm(canrange2,mean=priorR_mn2,sd=priorR_sd2,log=TRUE)
      
      MH2 <- canll2-curll2+prior_canll2-prior_curll2#+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
      
      if (log(runif(1))<MH2)
      {
        lrangeV=canrange2
        arv[ttt]=arv[ttt]+1
      }
      rangeU = exp(lrangeU)
      rangeV = exp(lrangeV)
      
      
      
      ##############################################:
      #####        KEEP TRACK OF STUFF       #######:
      ##############################################:
      
      keep.rangeU[iter,ttt] = rangeU
      keep.rangeV[iter,ttt] = rangeV
      keep.sigmaU[, iter,ttt] = sigmaU
      keep.sigmaV[, iter,ttt] = sigmaV
      # keep.taus1[iter,ttt] = taus1
      # keep.taus2[iter,ttt] = taus2
      keep.taue1[iter,ttt] = taue1
      keep.taue2[iter,ttt] = taue2
      keep.A[, iter,ttt] = A
      
      # keep.u1[,,iter,ttt]=as.matrix(U1)
      # keep.u2[,,iter,ttt]=as.matrix(U2)
      # keep.v2[,,iter,ttt]=as.matrix(V2)
      
      #keep.Y1.M[,,iter,ttt]=as.matrix(Y1)
      #keep.Y2.M[,,iter,ttt]=as.matrix(Y2)
      
      ##############################################:
      #####           PREDICTIONS            #######:
      ##############################################:
      
      if(iter>burn & predictions)# & iter %% 10 == 0
      {
        
        Mp=exp_corr(all.d, range=rangeU)
        Mp00=Mp[1:(n1+n2),1:(n1+n2)]
        Mp11=Mp[(n1+n2+1):(n1+n2+np1),(n1+n2+1):(n1+n2+np1)]
        Mp10=Mp[(n1+n2+1):(n1+n2+np1),1:(n1+n2)]
        Mp01=t(Mp10)
        
        E00=eigen(Mp00)
        E00.G=E00$vectors
        E00.D=E00$values
        
        Mp00.inv=E00.G%*%diag(1/E00.D)%*%t(E00.G)
        
        AA=Mp10%*%Mp00.inv
        a=Mp10%*%Mp00.inv%*%t(Mp10)
        a=round(a,digits=7) # to avoid numercial underflow
        B=Mp11-a
        
        Uls=rbind(U1,U2)  
        for (r in 1:nt)
        {
          Au=AA%*%Uls[,r] 
          sigmaB=sigmaU[r]*B
          #print(sigmaB)
          Ul.pred=rmvnorm(1,mean=Au,sigma=sigmaB)
          
          U1p[,r]=Ul.pred[(1:np1)]
        }
        
        
        for(i in 1:np1){Z1p[i,] <- fft_real(U1p[i,],inverse=TRUE)}
        Y1.pred <- beta1+Z1p+rnorm(n=nt,sd=sqrt(taue1)) #beta1?
        keep.Y1.P[,,iter,ttt]=as.matrix(Y1.pred)
      
        # Save Mean and Variance
        E.u1=E.u1+U1/(iters-burn)
        E.u2=E.u2+U2/(iters-burn)
        E.v=E.v+V2/(iters-burn)
        
        E2.u1=E2.u1+U1*U1/(iters-burn)
        E2.u2=E2.u2+U2*U2/(iters-burn)
        E2.v=E2.v+V2*V2/(iters-burn)
        
        # corr Y1 and Y2
        cory1y2=A*sigmaU/sqrt((sigmaU+taus1)*(A*sigmaU+sigmaV+taus2))
        E.corrY1Y2=E.corrY1Y2+cory1y2/(iters-burn)
        E2.corrY1Y2=E2.corrY1Y2+cory1y2*cory1y2/(iters-burn)
        
      }
    } # End of thin
  }
  Mu1=E.u1
  Mu2=E.u2
  Mv=E.v
  
  Vu1=E2.u1-E.u1*E.u1
  Vu2=E2.u2-E.u2*E.u2
  Vv2=E2.v-E.v*E.v
  
  McorY1Y2=E.corrY1Y2
  VcorY1Y2=E2.corrY1Y2-E.corrY1Y2*E.corrY1Y2
  
  keep.Y1.P=keep.Y1.P[,,(burn:iters),]
  out=list(keep.rangeU,keep.rangeV,keep.sigmaU,keep.sigmaV,keep.taue1,keep.taue2,keep.A,
           aru,arv,keep.Y1.P,Mu1,Mu2,Mv,Vu1,Vu2,Vv2,McorY1Y2,VcorY1Y2)
  names(out)=c('rangeU','rangeV','sigmaU','sigmaV','tau1','tau2','A','aru','arv','Y1.p','MeanU1',
               'MeanU2','MeanV2','VarU1','VarU2','VarV2','MeanCorrY1Y2','VarCorrY1Y2')
  return(out)
}
