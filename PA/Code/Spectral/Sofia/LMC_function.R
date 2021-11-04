
LMC_fit=function(Y1,Y2, s1,s2,
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
# keep.u1= array(0,dim=c(n1,nt,iters))
# keep.u2= array(0,dim=c(n2,nt,iters))
# keep.v2= array(0,dim=c(n2,nt,iters))
keep.Y1.M= array(0,dim=c(n1,nt,iters))
keep.Y2.M= array(0,dim=c(n2,nt,iters))

# start MCMC
start = proc.time()[3]

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
    # Ru should be a matrix and use data from all spectrum
    # Sweep operation to get rid of variance
    Ru=sweep(rbind(U1,U2),2,FUN='/',sqrt(sigmaU))
    Rv=sweep(V2,2,FUN='/',sqrt(sigmaV))
    
    #Ru = as.vector(U_sim)/sqrt(sigmaU[r])
    #Rv = as.vector(V2[,r])/sqrt(sigmaV[r])
    
    # range1
    Ms=exp_corr(d,range=exp(lrangeU))
    curll = sum(apply(Ru,2,dmvnorm,mean=rep(0,n1+n2),sigma=Ms,log=TRUE))
    #curll = sum(dmvnorm(t(Ru), mean=rep(0,n1+n2),sigma=Ms,log=TRUE))
    #curll = dmvnorm(Ru,rep(0,n1+n2),Ms,log=TRUE)
    canrange1 = rnorm(1,lrangeU,0.5)
    canM = exp_corr(d,range=exp(canrange1))
    canll = sum(apply(Ru,2,dmvnorm,mean=rep(0,n1+n2),sigma=canM,log=TRUE))
    #canll = sum(dmvnorm(t(Ru), mean=rep(0,n1+n2),sigma=canM,log=TRUE))
    #canll = dmvnorm(Ru,rep(0,n1+n2),canM,log=TRUE)
    
    MH1 <- canll-curll+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
    
    if (log(runif(1))<MH1)
    {
      lrangeU=canrange1
    }
    
    # range2
    Ss=exp_corr(dv2,range = exp(lrangeV))
    #curll2 = dmvnorm(Rv,rep(0,n2),Ss,log=TRUE)
    #curll2 = sum(apply(Rv,2,dmvnorm,mean=rep(0,n2),sigma=Ss,log=TRUE))
    curll2 = sum(dmvnorm(t(Rv), mean = rep(0,n2),sigma=Ss,log=TRUE ))
    canrange2 = rnorm(1,lrangeV,0.5)
    canS = exp_corr(dv2,range=exp(canrange2))
    #canll2 = dmvnorm(Rv,rep(0,n2),canS,log=TRUE)
    canll2 = sum(dmvnorm(t(Rv), mean = rep(0,n2),sigma=canS,log=TRUE ))
    
    MH2 <- canll2-curll2+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
    
    if (log(runif(1))<MH2)
    {
      lrangeV=canrange2
    }
    rangeU = exp(lrangeU)
    rangeV = exp(lrangeV)
    } # Ene of thin
  # Update range parameter

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
  
  # keep.u1[,,iter]=U1
  # keep.u2[,,iter]=U2
  # keep.v2[,,iter]=V2
  
  keep.Y1.M[,,iter]=as.matrix(Y1)
  keep.Y2.M[,,iter]=as.matrix(Y2)
}
proc.time()[3] - start


out=list(keep.rangeU,keep.rangeV,keep.sigmaU,keep.sigmaV,keep.taue1,keep.taue2,keep.A,keep.Y1.M,keep.Y2.M)
return(out)
}