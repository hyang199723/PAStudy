# Given longitude and latitude, output second-order location matrix 
# including intercept
locOpe <- function(lon, lat) {
  lon2 <- lon * lon
  lat2 <- lat * lat
  lonlat <- lon * lat
  o <- cbind(1, lon, lat, lon2, lat2, lonlat)
  return(o)
}

# Normalize one column
colNorm <- function(col) {
  return ((col - mean(col, na.rm = TRUE)) / sd(col, na.rm = TRUE))
}
# Normalize a given data frame columnwise
dfColNorm <- function(df) {
  df <- apply(df, 2, colNorm)
  return(df)
}

# Get start and end time for each month in a calendar year
getMonthTime <- function() {
  jan.start <- as.POSIXct('2020-01-01 00:00:00')
  feb.start <- as.POSIXct('2020-02-01 00:00:00')
  mar.start <- as.POSIXct('2020-03-01 00:00:00')
  apr.start <- as.POSIXct('2020-04-01 00:00:00')
  may.start <- as.POSIXct('2020-05-01 00:00:00')
  jun.start <- as.POSIXct('2020-06-01 00:00:00')
  jul.start <- as.POSIXct('2020-07-01 00:00:00')
  aug.start <- as.POSIXct('2020-08-01 00:00:00')
  sep.start <- as.POSIXct('2020-09-01 00:00:00')
  oct.start <- as.POSIXct('2020-10-01 00:00:00')
  nov.start <- as.POSIXct('2020-11-01 00:00:00')
  dec.start <- as.POSIXct('2020-12-01 00:00:00')
  
  cm.iterator <- c(jan.start, feb.start, mar.start, apr.start, may.start, jun.start,
                   jul.start, aug.start, sep.start, oct.start, nov.start, dec.start)
  
  return(cm.iterator)
}

# Read both 2020 EPA and PA data from Formatted_PA_FRM folder
readData <- function() {
  library(lubridate)
  pa <- read.csv('Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv')
  epa <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv')
  pa$Timestamp <- force_tz(as.POSIXct(pa$Timestamp), 'UTC')
  epa$Timestamp <- force_tz(as.POSIXct(epa$Timestamp), 'UTC')
  return(pa, epa)
}


cov.LMC <- function(d,t,theta){
  # t[i] = 1 means AQS; t[i] = 2 mean PA
  rho    <- theta[1]
  range1 <- exp(theta[2])
  range2 <- exp(theta[3])
  tau1   <- exp(theta[4])
  tau2   <- exp(theta[5])
  sig1   <- exp(theta[6])
  sig2   <- exp(theta[7])
  S      <- diag(ifelse(t==1,tau1,tau2)[,1])
  S[t==1,t==1] <- S[t==1,t==1] + sig1*exp(-d[t==1,t==1]/range1) +   
    sig2*rho*rho*exp(-d[t==1,t==1]/range2)
  S[t==1,t==2] <- S[t==1,t==2] + sig2*rho*exp(-d[t==1,t==2]/range2)
  S[t==2,t==1] <- S[t==2,t==1] + sig2*rho*exp(-d[t==2,t==1]/range2)
  S[t==2,t==2] <- S[t==2,t==2] + sig2*exp(-d[t==2,t==2]/range2)
  
  Sinv   <- solve(S)
  logdet <- determinant(S)$modulus[1]
  out    <- list(S=S,Sinv=Sinv,logdet=logdet)
  return(out)}   

# The log likelihood function
log_like <- function(Y,SIG){
  -0.5*SIG$logdet - 0.5*sum(Y*(SIG$Sinv%*%Y))
}

# The main MCMC function
LMC <- function(Y,X,s,type,
                mean_beta=0,sd_beta=10,
                mean_range=0,sd_range=1,
                mean_var=0,sd_var=1,
                mean_rho=0,sd_rho=10,
                iters=5000,burn=1000,thin=1,update=10){
  
  library(fields) 
  library(geoR)
  
  tick <- proc.time()[3]
  
  # Bookkeeping
  n        <- length(Y)
  p        <- ncol(X)
  
  d        <- as.matrix(dist(s))
  theta.mn <- c(mean_rho,mean_range,mean_range,
                mean_var,mean_var,mean_var,mean_var)
  theta.sd <- c(sd_rho,sd_range,sd_range,
                sd_var,sd_var,sd_var,sd_var)
  q        <- length(theta.mn)
  
  # Initial values
  
  theta <- theta.mn
  beta  <- lm(Y~X-1)$coef
  Xbeta <- X%*%beta
  Sig   <- cov.LMC(d,type,theta)
  curll <- log_like(Y-Xbeta,Sig)
  
  # Keep track of stuff
  
  keep_beta   <- matrix(0,iters,p)
  keep_theta  <- matrix(0,iters,q)
  colnames(keep_beta)  <- colnames(X)
  colnames(keep_theta) <- c("rho","log range1","log range2","log tau1","log tau2","log sig1","log sig2")
  att <- acc <- MH <- rep(0.5,q)
  
  # GO!!!
  
  for(iter in 1:iters){
    
    for(ttt in 1:thin){
      
      ##############################################:
      ####         MEAN PARAMETERS (GIBS)      #####:
      ##############################################:
      
      VVV   <- solve(t(X)%*%Sig$Sinv%*%X + diag(p)/sd_beta^2)
      MMM   <- t(X)%*%Sig$Sinv%*%Y + mean_beta/(sd_beta^2)
      beta  <- as.vector(VVV%*%MMM + t(chol(VVV))%*%rnorm(p)) 
      Xbeta <- X%*%beta
      curll <- log_like(Y-Xbeta,Sig)
      
      ##############################################:
      ####  COVARIANCE PARAMETERS (Metropolis) #####:
      ##############################################:
      
      for(j in 1:q){
        att[j]  <- att[j]+1
        cant    <- theta
        cant[j] <- rnorm(1,theta[j],MH[j])
        canSig  <- cov.LMC(d,type,cant)
        canll   <- log_like(Y-Xbeta,canSig)
        R       <- canll-curll+
          dnorm( cant[j],theta.mn[j],theta.sd[j],log=TRUE)-
          dnorm(theta[j],theta.mn[j],theta.sd[j],log=TRUE)
        if(log(runif(1))<R){
          acc[j] <- acc[j] + 1
          theta  <- cant
          Sig    <- canSig
          curll  <- canll
        }  
      }
      
      # Tuning
      for(j in 1:q){if(att[j]>50 & burn<iter){
        if(acc[j]/att[j]<0.2){MH[j] <- MH[j]*0.8}
        if(acc[j]/att[j]>0.6){MH[j] <- MH[j]*1.2}
        acc[j] <- att[j] <- 0
      }}
      
    } # end thin
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    
    keep_beta[iter,]  <- beta
    keep_theta[iter,] <- theta
    
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



# Do prediction at locations.
# theta: fitted theta_hat from LMC
# s: training locations
# s0: prediction locations
# x: training covaraites
# x0: prediction covariates
# type: training type (this should be the same as cov.LMC type)
# beta: fitted beta_hat from LMC
# y: True values at training locations
# Return: n * 1 predictted value, where n = nrow(s0)
LMCpredict <- function(theta, s, s0, x, x0, type, beta, y) {
  d <- as.matrix(dist(s))
  sigma1_inv <- cov.LMC(d, type, theta)$Sinv
  
  # Construct sigma0
  rho    <- theta[1]
  range1 <- exp(theta[2])
  range2 <- exp(theta[3])
  tau1   <- exp(theta[4])
  tau2   <- exp(theta[5])
  sig1   <- exp(theta[6])
  sig2   <- exp(theta[7])
  
  s1 <- rbind(s0, s)
  # distance matrix between s0 and s
  dist <- as.matrix(dist(s1))
  rownames(dist) <- NULL
  
  d <- as.matrix(dist[1:nrow(s0), (nrow(s0)+1):ncol(dist)])
  if (nrow(s0) == 1) {d <- t(d)}
  
  S      <- matrix(0,nrow = nrow(s0), ncol = ncol(d))
  S[,type==1] <- sig1*exp(-d[,type==1]/range1) +   
    sig2*rho*rho*exp(-d[,type==1]/range2)
  S[,type==2] <- sig2*rho*exp(-d[,type==2]/range2)
  sigma0 <- S
  
  y_hat <- x0%*%beta + sigma0 %*% sigma1_inv %*% (y - x %*% beta)
  return(y_hat)
}


LMCcovariance <- function(theta, s, s0, x, x0, type, beta, y) {
  d <- as.matrix(dist(s))
  sigma1_inv <- cov.LMC(d, type, theta)$Sinv
  
  # Construct sigma0
  rho    <- theta[1]
  range1 <- exp(theta[2])
  range2 <- exp(theta[3])
  tau1   <- exp(theta[4])
  tau2   <- exp(theta[5])
  sig1   <- exp(theta[6])
  sig2   <- exp(theta[7])
  
  s1 <- rbind(s0, s)
  # distance matrix between s0 and s
  dist <- as.matrix(dist(s1))
  rownames(dist) <- NULL
  
  d <- as.matrix(dist[1:nrow(s0), (nrow(s0)+1):ncol(dist)])
  if (nrow(s0) == 1) {d <- t(d)}
  
  S      <- matrix(0,nrow = nrow(s0), ncol = ncol(d))
  S[,type==1] <- sig1*exp(-d[,type==1]/range1) +   
    sig2*rho*rho*exp(-d[,type==1]/range2)
  S[,type==2] <- sig2*rho*exp(-d[,type==2]/range2)
  sigma0 <- S
  var1 <- sig1 + (rho^2) * sig2 + tau1
  var <- matrix(var1, nrow = nrow(s0), ncol = 1) -
    diag(diag(sigma0 %*% sigma1_inv %*% t(sigma0)), nrow = nrow(s0), ncol = 1)
  return(var)
}





































