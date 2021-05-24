#######################################################
#
# INPUTS
#
# y          := n-1 of observations
# s          := nx2 matrix of spatial coordinates
# X          := nxp matrix of covariates
#
# sp,Xp      := s and X for prediction locations
#
# mean_r     := r = logit(parsill/(var)) ~ N(mean_r,sd_r)
# sd_r
# mean_nu    := log(nu) ~ N(mean_nu,sd_nu)
# sd_nu
# mean_range := log(rangee) ~ N(mean_range,sd_range)
# sd_range
# a_var      := variance ~ InvGamma(a_var,b_var)
# sd_beta    := beta[j] ~ N(0,sd_beta)
#
# iters      := number of MCMC iterations
# burn       := length of burn-in
#
#######################################################



Bayes_Krige<-function(y,s,X,
                      sp=NULL,Xp=NULL,
                      mean_r=0,sd_r=1,
                      mean_nu=-1,sd_nu=1,
                      mean_range=0,sd_range=10,
                      a_var=.01,b_var=.01,
                      sd_beta=1000,
                      iters=5000,burn=1000){
  
  library(mvtnorm)
  
  # Bookkeeping
  
  n        <- length(y)
  p        <- ncol(X)
  d        <- rdist(s,s)
  diag(d)  <- 0
  theta.mn <- c(mean_r,mean_range,mean_nu)
  theta.sd <- c(sd_r,sd_range,sd_nu)
  
  predictions<-!is.null(sp) & !is.null(Xp)
  np<-1
  if(predictions){
    np        <- nrow(sp)
    d12       <- rdist(sp,s)
    d11       <- rdist(sp,sp)
    diag(d11) <- 0
  }
  
  # Initial values
  
  beta    <- lm(y~X-1)$coef
  sigma2  <- var(as.vector(y-X%*%beta))
  theta   <- c(0,0,0)   
  Sigma   <- corfx(d,theta)
  Siginv  <- solve(Sigma)
  
  # Keep track of stuff
  
  keep.beta <- matrix(0,iters,p)
  keepers   <- matrix(0,iters,4)
  pred      <- matrix(0,iters,np)
  
  colnames(keepers)   <- c("sigma2","r","range","nu")
  colnames(keep.beta) <- colnames(X)
  
  # GO!!!
  
  for(i in 1:iters){
    
    ##############################################:
    #####       MEAN PARAMETERS (Gibbs)    #######:
    ##############################################:
    
    tXS  <- t(X)%*%Siginv/sigma2    
    VVV  <- solve(tXS%*%X + diag(p)/sd_beta^2)
    MMM  <- VVV%*%tXS%*%y
    beta <- MMM + t(chol(VVV))%*%rnorm(p)
    
    ##############################################:
    #####          VARIANCE (Gibbs)        #######:
    ##############################################:
    
    R      <- y-X%*%beta
    a      <- n/2+a_var
    b      <- t(R)%*%Siginv%*%R/2+b_var
    sigma2 <- 1/rgamma(1,a,b) 
    
    ##############################################:
    #### CORRELATION PARAMETERS (Metropolis) #####:
    ##############################################:
    
    R      <- as.vector(y-X%*%beta)/sqrt(sigma2)
    curll  <- dmvnorm(R,rep(0,n),Sigma,log=TRUE)
    
    for(j in 1:3){
      cantheta    <- theta
      cantheta[j] <- rnorm(1,theta[j],0.5)
      canSigma    <- corfx(d,cantheta)
      canll       <- dmvnorm(R,rep(0,n),canSigma,log=TRUE)
      
      MH <- canll-curll+
        dnorm(cantheta[j],theta.mn[j],theta.sd[j],log=TRUE)-
        dnorm(   theta[j],theta.mn[j],theta.sd[j],log=TRUE)
      if(log(runif(1))<MH){
        theta  <- cantheta
        Sigma  <- canSigma
        curll  <- canll
      }  
    }
    Siginv <- solve(Sigma)
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep.beta[i,]  <- beta
    keepers[i,1]   <- sigma2
    keepers[i,2]   <- 1/(1+exp(-theta[1]))
    keepers[i,3:4] <- exp(theta[2:3])
    
    ##############################################:
    #####           PREDICTIONS            #######:
    ##############################################:
    
    if(i>burn & predictions){
      S11      <- corfx(d11,theta)*sigma2
      S12      <- corfx(d12,theta)*sigma2
      S22inv   <- Siginv/sigma2
      VVV      <- S11-S12%*%S22inv%*%t(S12)
      MMM      <- Xp%*%beta + 
        S12%*%S22inv%*%(y-X%*%beta)
      pred[i,] <- MMM + t(chol(VVV))%*%rnorm(np)
    }
    
  }   
  
  
  list(beta=keep.beta,keepers=keepers,pred=pred)}


corfx <- function(d,theta){
  r   <- 1/(1+exp(-theta[1])) # sigma-to-noise ratio
  rho <- exp(theta[2])        # range
  nu  <- exp(theta[3])        # smoothness
  COR <- r*matern(d,rho,nu)
  COR[d==0]<-1
  return(COR)}


### Simulate a fake dataset
library(fields)
library(geoR)

n     <- 200
p     <- 4
sig2  <- 3
r     <- 0.8
rho   <- 0.1
nu    <- 0.5
beta  <- c(10,0,1,2)
theta <- c(log(r/(1-r)),log(rho),log(nu))   

S     <- cbind(runif(n),runif(n))
X     <- matrix(rnorm(n*p),n,p)
X[,1] <- 1   
d     <- rdist(S);diag(d)<-0
Cor   <- corfx(d,theta)
Y     <- X%*%beta + t(chol(Cor))%*%rnorm(n,0,sqrt(sig2))

test <- runif(n)<0.1
S1   <- S[!test,]
X1   <- X[!test,]
Y1   <- Y[!test]
S2   <- S[test,]
X2   <- X[test,]
Y2   <- Y[test]


fit <- Bayes_Krige(Y1,S1,X1,sp=S2,Xp=X2)


hist(fit$keepers[1000:5000,1],main="Variance",breaks=25)
abline(v=sig2,lwd=2,col=4)

hist(fit$keepers[1000:5000,2],main="Signal-to-noise",breaks=25)
abline(v=r,lwd=2,col=4)

hist(fit$keepers[1000:5000,3],main="Range",breaks=25)
abline(v=rho,lwd=2,col=4)

hist(fit$keepers[1000:5000,4],main="Smoothness",breaks=25)
abline(v=nu,lwd=2,col=4)

pairs(fit$keepers[1000:5000,])

for(j in 1:p){
  hist(fit$beta[1000:5000,j],main=paste0("beta",j),breaks=25)
  abline(v=beta[j],lwd=2,col=4)
}

boxplot(fit$pred[1000:5000,],main="Predictions",outline=FALSE)
lines(Y2,col=4,lwd=2)