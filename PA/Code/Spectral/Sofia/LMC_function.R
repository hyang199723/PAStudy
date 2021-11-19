# This includes the main LMC functions
# This script has two functions, LMC_fit and compact.LMC_fit
# Last update: 11/18/2021
LMC_fit=function(Y1,Y2, s1,s2,
                 mean_range=0, sd_range=1, mean_var=0, sd_var=1, mean_rho=0,
                 sd_rho=10, iters=3000, burn=1000, thin=1, update=10) {
  n1       <- nrow(Y1)
  n2       <- nrow(Y2)
  nt       <- ncol(Y1)
  m1       <- is.na(Y1)
  m2       <- is.na(Y2)
  d        <- as.matrix(dist(rbind(s1,s2)))
  dv2 = as.matrix(dist(s2))
  const    <- nt/2
  
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
  keep.taus1 = matrix(0, iters,thin)
  keep.taus2 = matrix(0, iters,thin)
  keep.taue1 = matrix(0, iters,thin)
  keep.taue2 = matrix(0, iters,thin)
  # Vector parameters; Dimension of data structure: #stations * nt * iters
  # keep.u1= array(0,dim=c(n1,nt,iters))
  # keep.u2= array(0,dim=c(n2,nt,iters))
  # keep.v2= array(0,dim=c(n2,nt,iters))
  keep.Y1.M= array(0,dim=c(n1,nt,iters,thin))
  keep.Y2.M= array(0,dim=c(n2,nt,iters,thin))


  ## get info for the ranges priors 

  priorR_mn1 <- log(max(d)) - 1.5
  priorR_sd1 <- 1

  priorR_mn2 <- log(max(dv2)) - 1.5
  priorR_sd2 <- 1



  # start MCMC
  for(ttt in 1:thin){
    # Mean imputation
    beta1 <- mean(colMeans(Y1,na.rm=TRUE))
    beta2 <- mean(colMeans(Y2,na.rm=TRUE))
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
    rangeU  <- 0.2
    rangeV   <- 0.2
    lrangeU  <- log(rangeU)
    lrangeV   <-log(rangeV)
    taue1  <- var(as.vector(Y1),na.rm=TRUE)
    taue2  <- var(as.vector(Y2),na.rm=TRUE)
    sigmaU   <- rep(taue1,nt)
    sigmaV   <- rep(taue1,nt)
    A      <- rep(0,nt)
    
    
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
      
      # Gibbs sampling
      for (r in 1:nt) # for each spectral 
      {
        # # Sample U1
        newDiag = 1/taus1 + 1/sigmaU[r] * 1/S1_D
        sigmaU1 <- S1_G %*% diag(1/newDiag) %*% t(S1_G)
        # S1_G
        meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1invA1 %*% U2[,r])
        U1[,r] = meanU1+S1_G%*%rnorm(n1, 0, 1/sqrt(newDiag))
        
        # # Sample U2
        newDiag = A[r]^2/taus2 + 1/sigmaU[r] * 1/S2_D
        sigmaU2 <- S2_G %*% diag(1/newDiag) %*% t(S2_G)
        
        meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 *A[r] * V2[,r] +
                                 1/sigmaU[r] * S2invA2 %*% U1[,r])
        U2[,r] = meanU2+S2_G%*%rnorm(n2,0,1/sqrt(newDiag))
        
        # # Sample V2
        newDiag = 1/sigmaV[r] * 1/V_D + 1/taus2
        sigmaV2 = V_G %*% diag(1/newDiag) %*% t(V_G)
        
        meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
        V2[,r] =meanV2+V_G%*%rnorm(n2,0,1/sqrt(newDiag))      
        
        # # Sample Al
        sigmaAl <- solve(1/taus2 * t(U2[,r]) %*% U2[,r] + 1/5)
        meanAl <- sigmaAl %*% (1/taus2 * t(Ys2[,r]) %*% U2[,r] - 1/taus2 * t(V2[,r]) %*% U2[,r] + 0.8/5)
        A[r] <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl, sd=sqrt(sigmaAl))
        
        # # Sample sig1
        U_sim <- as.vector(append(U1[,r], U2[,r]))
        a <- (n2+n1)/2 + 1
        b <- (t(U_sim) %*% eigU$Q %*% U_sim) / 2 + 1
        sigmaU[r] = 1/rgamma(1, a, b)
        
        # # # Sample sig2
        a <- n2/2 + 1
        b <- (t(V2[,r]) %*% eigV$Q %*% V2[,r]) / 2 + 1 
        sigmaV[r] = 1/rgamma(1, a, b)
        
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
      
      MH1 <- canll-curll+prior_canll-prior_curll+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
      
      if (log(runif(1))<MH1)
      {
        lrangeU=canrange1
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
      
      MH2 <- canll2-curll2+prior_canll2-prior_curll2+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
      
      if (log(runif(1))<MH2)
      {
        lrangeV=canrange2
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
      keep.taus1[iter,ttt] = taus1
      keep.taus2[iter,ttt] = taus2
      keep.taue1[iter,ttt] = taue1
      keep.taue2[iter,ttt] = taue2
      keep.A[, iter,ttt] = A
      
      # keep.u1[,,iter]=U1
      # keep.u2[,,iter]=U2
      # keep.v2[,,iter]=V2
      
      keep.Y1.M[,,iter,ttt]=as.matrix(Y1)
      keep.Y2.M[,,iter,ttt]=as.matrix(Y2)
      #print(iter)
    } # End of thin
  }
  out=list(keep.rangeU,keep.rangeV,keep.sigmaU,keep.sigmaV,keep.taue1,keep.taue2,keep.A,keep.Y1.M,keep.Y2.M)
  names(out)=c('rangeU','rangeV','sigmaU','sigmaV','tau1','tau2','A','Y1.m','Y2.m')
  return(out)
}



compact.LMC_fit=function(Y1,Y2, s1,s2,
                 mean_range=0, sd_range=1, mean_var=0, sd_var=1, mean_rho=0,
                 sd_rho=10, iters=3000, burn=1000, thin=1, update=10) {
  n1       <- nrow(Y1)
  n2       <- nrow(Y2)
  nt       <- ncol(Y1)
  m1       <- is.na(Y1)
  m2       <- is.na(Y2)
  d        <- as.matrix(dist(rbind(s1,s2)))
  dv2 = as.matrix(dist(s2))
  const    <- nt/2
  
  # Mean imputation
  beta1 <- mean(colMeans(Y1,na.rm=TRUE))
  beta2 <- mean(colMeans(Y2,na.rm=TRUE))
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
  rangeU  <- 0.2
  rangeV   <- 0.2
  lrangeU  <- log(rangeU)
  lrangeV   <-log(rangeV)
  taue1  <- var(as.vector(Y1),na.rm=TRUE)
  taue2  <- var(as.vector(Y2),na.rm=TRUE)
  sigmaU   <- rep(taue1,nt)
  sigmaV   <- rep(taue1,nt)
  A      <- rep(0,nt)
  
  # Keep track of stuff
  
  # All parameters: "rangeu","rangev","sigmaU","sigmaV","taue1","taue2","A",'beta1','beta2'
  
  # Parameters that are different for each spectrum: sigmaU, sigmaV, A, beta1, beta2
  # Dimension of data structure: nt * iterations
  keep.sigmaU = array(0, dim = c(nt, iters))
  keep.sigmaV = array(0, dim = c(nt, iters))
  keep.A = array(0, dim = c(nt, iters))
  
  # Parameters that stay the same: rangeU, rangeV, taus1, taus2
  # Dimension of data structure: 1 * iterations
  keep.rangeU = rep(0, iters)
  keep.rangeV = rep(0, iters)
  keep.taus1 = rep(0, iters)
  keep.taus2 = rep(0, iters)
  keep.taue1 = rep(0, iters)
  keep.taue2 = rep(0, iters)
  # Vector parameters; Dimension of data structure: #stations * nt * iters
  keep.u1= array(0,dim=c(n1,nt,iters))
  keep.u2= array(0,dim=c(n2,nt,iters))
  keep.v2= array(0,dim=c(n2,nt,iters))
  keep.Y1.M= array(0,dim=c(n1,nt,iters))
  keep.Y2.M= array(0,dim=c(n2,nt,iters))
  
  ## get info for the ranges priors 
  
  priorR_mn1 <- log(max(d)) - 1.5
  priorR_sd1 <- 1
  
  priorR_mn2 <- log(max(dv2)) - 1.5
  priorR_sd2 <- 1

  # set some constant
  a1 <- (n2+n1)/2 + 1
  a2 <- n2/2 + 1
  
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
      
      # Gibbs sampling
      
      E = matrix(rnorm(n1*n2),n1,nt)
      E2 = matrix(rnorm(n1*n2),n2,nt)
      E3 = matrix(rnorm(n1*n2),n2,nt)
      
      
      GYs1 <- t(S1_G)%*%Ys1
      GVYs2 <- t(V_G)%*%Ys2
      GYs2 <- t(S2_G)%*%Ys2
      
      GSU2 <- t(S1_G)%*%S1invA1%*%U2
      
      # Sample U1
      ND1=t(1/taus1 +apply(as.matrix(1/S1_D),1,FUN='*',1/sigmaU))
      ## Mean part with sweep function 
      Mfs=GYs1/taus1+sweep(GSU2,2,FUN='/',sigmaU)
      Mean.P=sweep(Mfs,1,FUN='/',ND1)
      U1 = S1_G%*%(Mean.P+E)
      
      
      GSU1 <- t(S2_G)%*%S2invA2%*%U1
      GSV2 <- t(S2_G)%*%V2
      
      
      # Sample U2
      ND2=t(A^2/taus2 +apply(as.matrix(1/S2_D),1,FUN='*',1/sigmaU))
      Mfs2=sweep((GYs2),1,FUN='*',A)/taus2- sweep(GSV2,2,FUN='*',A)/taus2+sweep(GSU1,2,FUN='/',sigmaU)
      Mean.P2=sweep(Mfs2,1,FUN='/',ND2)
      U2 = S2_G%*%(Mean.P2+E2)

      GVU2 <- t(V_G)%*%U2
      U_sim=rbind(U1,U2)

      # Sample V2
      ND3=t(apply(as.matrix(1/V_D),1,FUN='*',1/sigmaV))+1/taus2
      Mfs3=GVYs2/taus2- sweep(GVU2,1,FUN='*',A)/taus2
      Mean.P3=sweep(Mfs3,1,FUN='/',ND3)
      V2 = V_G%*%(Mean.P3+E3)

      # Sample A
      U2U2.taus=diag(t(U2)%*%U2)/taus2+1/5
      sigmaAl=1/(U2U2.taus)
      p2=diag(t(Ys2)%*%U2)/taus2-diag(t(V2)%*%U2)/taus2+0.8/5
      meanAl=sigmaAl*p2
      for (r in 1:nt){A[r] <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl[r], sd=sqrt(sigmaAl[r]))}
      
      #Sample sigmaU1
      b1 <- (t(U_sim) %*% eigU$Q %*% U_sim) / 2 + 1
      for(r in 1:nt){sigmaU[r] = 1/rgamma(1, a1, diag(b1)[r])}
      
      #Sample sigmaU2
      b2<- (t(V2) %*% eigV$Q %*% V2) / 2 + 1 
      for(r in 1:nt){sigmaV[r] = 1/rgamma(1, a2, diag(b2)[r])}
      
      #for (r in 1:nt) # for each spectral 
      #{
        # # Sample U1
        # newDiag = 1/taus1 + 1/sigmaU[r] * 1/S1_D
        # sigmaU1 <- S1_G %*% diag(1/newDiag) %*% t(S1_G)
        # # S1_G
        # meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1invA1 %*% U2[,r])
        # U1[,r] = meanU1+S1_G%*%rnorm(n1, 0, 1/sqrt(newDiag))
        
        # # Sample U2
        # newDiag = A[r]^2/taus2 + 1/sigmaU[r] * 1/S2_D
        # sigmaU2 <- S2_G %*% diag(1/newDiag) %*% t(S2_G)
        # 
        # meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 *A[r] * V2[,r] +
        #                          1/sigmaU[r] * S2invA2 %*% U1[,r])
        # U2[,r] = meanU2+S2_G%*%rnorm(n2,0,1/sqrt(newDiag))

        # # Sample V2
        # newDiag = 1/sigmaV[r] * 1/V_D + 1/taus2
        # sigmaV2 = V_G %*% diag(1/newDiag) %*% t(V_G)
        # 
        # meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
        # V2[,r] =meanV2+V_G%*%rnorm(n2,0,1/sqrt(newDiag))

        # # Sample Al
        # sigmaAl <- solve(1/taus2 * t(U2[,r]) %*% U2[,r] + 1/5)
        # meanAl <- sigmaAl %*% (1/taus2 * t(Ys2[,r]) %*% U2[,r] - 1/taus2 * t(V2[,r]) %*% U2[,r] + 0.8/5)
        # A[r] <- rtruncnorm(1, a=0, b=+Inf, mean=meanAl, sd=sqrt(sigmaAl))
        
        # # Sample sig1
        #U_sim <- as.vector(append(U1[,r], U2[,r]))
        #b <- (t(U_sim[,r]) %*% eigU$Q %*% U_sim[,r]) / 2 + 1
        #sigmaU[r] = 1/rgamma(1, a1, b)
        
        
        
        
        # # # Sample sig2
        # b <- (t(V2[,r]) %*% eigV$Q %*% V2[,r]) / 2 + 1 
        # sigmaV[r] = 1/rgamma(1, a2, b)
      #}
      
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
      
      MH1 <- canll-curll+prior_canll-prior_curll+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
      
      if (log(runif(1))<MH1)
      {
        lrangeU=canrange1
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
      
      MH2 <- canll2-curll2+prior_canll2-prior_curll2+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
      
      if (log(runif(1))<MH2)
      {
        lrangeV=canrange2
      }
      rangeU = exp(lrangeU)
      rangeV = exp(lrangeV)
      
      
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep.rangeU[iter] = rangeU
    keep.rangeV[iter] = rangeV
    keep.sigmaU[, iter] = sigmaU
    keep.sigmaV[, iter] = sigmaV
    keep.taus1[iter] = taus1
    keep.taus2[iter] = taus2
    keep.taue1[iter] = taue1
    keep.taue2[iter] = taue2
    keep.A[, iter] = A
    
    keep.u1[,,iter]=U1
    keep.u2[,,iter]=U2
    keep.v2[,,iter]=V2

    keep.Y1.M[,,iter]=as.matrix(Y1)
    keep.Y2.M[,,iter]=as.matrix(Y2)
    #print(iter)
  }
  out=list(keep.rangeU,keep.rangeV,keep.sigmaU,keep.sigmaV,keep.taue1,keep.taue2,keep.A,keep.Y1.M,keep.Y2.M)
  names(out)=c('rangeU','rangeV','sigmaU','sigmaV','tau1','tau2','A','Y1.m','Y2.m')
  return(out)
  }
}