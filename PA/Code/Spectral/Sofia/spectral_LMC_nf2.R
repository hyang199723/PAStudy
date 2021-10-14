#######################################################
# same range for all sprectrums 
#######################################################

source('ExtraFunctions.R')
source('simAllTS.R')
library(fields) 
library(geoR)
library(truncnorm)

# set initial values 
mean_range=0
sd_range=1
mean_var=0
sd_var=1
mean_rho=0
sd_rho=10
iters=3000
burn=1000
thin=1
update=10

# Bookkeeping
n1       <- nrow(Y1)
n2       <- nrow(Y2)
nt       <- ncol(Y1)
m1       <- is.na(Y1)
m2       <- is.na(Y2)
d        <- as.matrix(dist(rbind(s1,s2)))
dv2=as.matrix(dist(coords[type==2,]))
const    <- nt/2# ?

# Initial values
beta1 <- mean(Y1,na.rm=TRUE)
beta2 <- mean(Y2,na.rm=TRUE)
Y1[m1] <- beta1
Y2[m2] <- beta2
U1     <- matrix(0,n1,nt)
U2     <- matrix(0,n2,nt)
V2     <- matrix(0,n2,nt)
rangeU  <- exp(mean_rho)
rangeV   <- exp(mean_rho)
lrangeU  <- log(rangeU)
lrangeV   <-log(rangeV)

taue1  <- var(as.vector(Y1))
taue2  <- var(as.vector(Y2))
sigmaU   <- rep(taue1,nt)
sigmaV   <- rep(taue1,nt)
A      <- rep(1,nt)
Z1     <- matrix(0,n1,nt)
Z2     <- matrix(0,n2,nt)
Ys1    <- matrix(0,n1,nt)
Ys2    <- matrix(0,n2,nt)

# Keep track of stuff

keep_theta  <- array(0,dim=c(iters,9,nt))
keep.u1= array(NA,dim=c(iters,n1,nt))
keep.u2= array(NA,dim=c(iters,n2,nt))
keep.v2= array(NA,dim=c(iters,n2,nt))
keep.Y1.M= array(NA,dim=c(iters,dim(m1)))
keep.Y2.M= array(NA,dim=c(iters,dim(m2)))

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
    
    Y1[m1] <- rnorm(sum(m1),beta1+Z1[m1],sqrt(taue1)) 
    Y2[m2] <- rnorm(sum(m2),beta2+Z2[m2],sqrt(taue2)) 
    
    ##############################################:
    ####      MEAN PARAMETERS (real space)   #####:
    ##############################################:
    
    # full conditional for beta1 and beta2 ...
    #VVV   <- (taue1*n1*nt + 0.01)
    #MMM   <- taue1*sum(Y1-Z1) 
    #beta1 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
    beta1=beta.1
    
    # VVV   <- (taue2*n2*nt + 0.01)
    # MMM   <- taue2*sum(Y2-Z2) 
    #beta2 <- rnorm(1,MMM/VVV,1/sqrt(VVV))
    beta2=beta.2
    
    ##############################################:
    ####     ERROR VARIANCES (real space)    #####:
    ##############################################:
    # full conditionals for taue1 and taue2 
    taue1 <- 1/rgamma(1,n1*nt/2+0.01,sum((Y1-beta1-Z1)^2)/2+0.01)
    taue2 <- 1/rgamma(1,n2*nt/2+0.01,sum((Y2-beta2-Z2)^2)/2+0.01)
    # 
    ##############################################:
    ####     Transform to spectral land      #####:
    ##############################################:
    
    for(i in 1:n1){Ys1[i,] <- fft_real(Y1[i,]-beta1)}
    for(i in 1:n2){Ys2[i,] <- fft_real(Y2[i,]-beta2)}
    taus1 <- const*taue1
    taus2 <- const*taue2
    
    ##############################################:
    ####      LMC TERMS (spectral space)     #####:
    ##############################################:
    
    eigU = invU(d,n1,n2,rangeU)
    S11 = eigU$S11; S12 <- eigU$S12; S21 <- eigU$S21; S22 <- eigU$S22
    S1 = S11 - eigU$A12 %*% S21
    ES1=eigen(S1)
    G   = ES1$vectors
    D   = ES1$values
    S1inv = G%*%diag(1/D)%*%t(G)
    A1 = S12 %*% solve(S22)
    
    S11inv=solve(S11)
    S2 = S22 - S21 %*% S11inv %*% S12
    A2 = S21 %*% S11inv
    ES2 = eigen(S2)
    G = ES2$vectors
    D = ES2$values
    S2inv = G%*%diag(1/D)%*%t(G)
    
    eigV  = invV(dv2,rangeV)
    
    for (r in 1:nt) # for each spectral 
    {
      # # Sample U1
      sigmaU1 <- solve(1/taus1 * diag(1, n1) + 1/sigmaU[r] * S1inv)
      meanU1 <- sigmaU1 %*% (1/taus1 * Ys1[,r] + 1/sigmaU[r] * S1inv %*% A1 %*% U2[,r])
      U1[,r] <- as.vector(t(chol(sigmaU1)) %*% rnorm(n1)) + meanU1
      
      # # Sample U2
      sigmaU2 <- solve(A[r]^2/taus2 * diag(1, n2) + 1/sigmaU[r] * S2inv)
      meanU2 <- sigmaU2 %*% (1/taus2 * A[r] * Ys2[,r] - 1/taus2 * A[r] * V2[,r] +
                               1/sigmaU[r] * S2inv %*% A2 %*% U1[,r])
      U2[,r] <- as.vector(t(chol(sigmaU2)) %*% rnorm(n2)) + meanU2
      
      # # Sample V2
      sigmaV2 <- solve(1/sigmaV[r] *eigV$Q + 1/taus2 * diag(1, n2))
      meanV2 <- sigmaV2 %*% (1/taus2 * Ys2[,r] - 1/taus2 * A[r] * U2[,r])
      V2[,r] <- as.vector(t(chol(sigmaV2)) %*% rnorm(n2)) + meanV2
      
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
      
      ###############################################
      ##       Metropolis H: Range parameters      ##
      ###############################################
      
      Ru = as.vector(U_sim)/sqrt(sigmaU[r])
      Rv = as.vector(V2[,r])/sqrt(sigmaV[r])
      
      # range1
      Ms=exp_corr(d,range=exp(lrangeU))
      curll = dmvnorm(Ru,rep(0,n1+n2),Ms,log=TRUE)
      canrange1 = rnorm(1,lrangeU,0.5)
      canM = exp_corr(d,range=exp(canrange1))
      canll = dmvnorm(Ru,rep(0,n1+n2),canM,log=TRUE)
      
      MH1 <- canll-curll+dnorm(canrange1,log=TRUE)-dnorm(lrangeU,log=TRUE)
      
      if (log(runif(1))<MH1)
      {
        lrangeU=canrange1
      }
      
      # range2
      Ss=exp_corr(dv2,range = exp(lrangeV))
      curll2 = dmvnorm(Rv,rep(0,n2),Ss,log=TRUE)
      canrange2 = rnorm(1,lrangeV,0.5)
      canS = exp_corr(dv2,range=exp(canrange2))
      canll2 = dmvnorm(Rv,rep(0,n2),canS,log=TRUE)
      
      MH2 <- canll2-curll2+dnorm(canrange2,log=TRUE)-dnorm(lrangeV,log=TRUE)
      
      if (log(runif(1))<MH2)
      {
        lrangeV=canrange2
      }
      
      
    }
    
    
    
  } # end thin
  
  ##############################################:
  #####        KEEP TRACK OF STUFF       #######:
  ##############################################:
  
  
  keep_theta[iter, 1,]  <- rep(exp(lrangeU),nt)
  keep_theta[iter, 2,]  <- rep(exp(lrangeV),nt)
  keep_theta[iter, 3,]  <- sigmaU
  keep_theta[iter, 4,]  <- sigmaV
  keep_theta[iter, 5,]  <- rep(taue1,nt) # because we have the same tau for all spectrums 
  keep_theta[iter, 6,]  <- rep(taue2,nt)
  keep_theta[iter, 7,]  <- A
  keep_theta[iter, 8,]  <- rep(beta1,nt)
  keep_theta[iter, 9,]  <- rep(beta2,nt)
  
  keep.u1[iter,,]=U1
  keep.u2[iter,,]=U2
  keep.v2[iter,,]=V2
  
  keep.Y1.M[iter,,]=Y1
  keep.Y2.M[iter,,]=Y2
  
  
}   


### Results 

# create data frame 
P=c("rangeu","rangev","sigmaU","sigmaV","taue1","taue2","A",'beta1','beta2')
rv=list(rep(rangeu,nt),rep(rangev,nt),sigmau,sigmav,rep(tau1,nt),rep(tau2,nt),al,
        rep(beta1,nt),rep(beta2,nt))
values=c()
param=c()
n.iter=c()
time=c()
realvalues=c()
for(i in 1:9)
{
  for (j in 1:nt)
  {
    values=c(values,keep_theta[,i,j])
    param=c(param,rep(P[i],iters))
    n.iter=c(n.iter,rep(1:iters))
    time=c(time,rep(j,iters))
    realvalues=c(realvalues,rep(rv[[i]][j],iters))
    
  }
}



mcmc.results=data.frame(n.iter,values,param,time,realvalues)
mcmc.results$n.iter=as.numeric(mcmc.results$n.iter)

#al 
ggplot(mcmc.results %>% filter(param=='A',time%in%seq(1:10)))+geom_line(aes(x=n.iter,y=values))+
  facet_wrap(~time, scales = "free")+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+
  ggtitle('Al values')
#tau1
ggplot(mcmc.results %>% filter(param=='taue1'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+ylim(c(0,5))+ggtitle('tau1 values')
#tau2
ggplot(mcmc.results %>% filter(param=='taue2'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+ggtitle('tau2 values')

#sigmau
ggplot(mcmc.results %>% filter(param=='sigmaU'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+ggtitle('sigma1 values')+
  facet_wrap(~time,scales = "free")

#sigmav
ggplot(mcmc.results %>% filter(param=='sigmaV'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+ggtitle('sigma1 values')+
  facet_wrap(~time,scales = "free")

#phou
ggplot(mcmc.results %>% filter(param=='rangeu'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+ggtitle('phou values')

ggplot(mcmc.results %>% filter(param=='rangev'))+geom_line(aes(x=n.iter,y=values))+theme_bw()+geom_hline(aes(yintercept = realvalues,col='red'))+
  ggtitle('phov')

# missing values 
burniterns=seq((iters/2),iters)
misY=c()
n.iter=c()
n.time=c()
n.site=c()
misstype=c()
n.time=c()
n.type=c()

for (t in 1:nt)
{
  misY=c(misY,c(keep.Y1.M[burniterns,,t]))
  n.site=c(n.site,rep(1:(dim(m1)[1]),each=length(burniterns)))
  n.iter=c(n.iter,rep(burniterns,dim(m1)[1]))
  misstype=c(misstype,rep(1,length(keep.Y1.M[burniterns,,t])))
  n.time=c(n.time,rep(t,length(c(keep.Y1.M[burniterns,,t]))))
  n.type=c(n.type,rep(1,length(c(keep.Y1.M[burniterns,,t]))))
  
  misY=c(misY,c(keep.Y2.M[burniterns,,t]))
  n.site=c(n.site,rep(1:(dim(m2)[1]),each=length(burniterns)))
  n.iter=c(n.iter,rep(burniterns,dim(m2)[1]))
  misstype=c(misstype,rep(1,length(keep.Y2.M[burniterns,,t])))
  n.time=c(n.time,rep(t,length(c(keep.Y2.M[burniterns,,t]))))
  n.type=c(n.type,rep(2,length(c(keep.Y2.M[burniterns,,t]))))
  
}


MY=data.frame(misY,n.site,misstype,n.iter,n.time,n.type)

#Calculate posterior Mean
pm=MY %>% group_by(n.site,n.time,n.type) %>% summarise(vals=mean(misY))
pm$TF=rep('estimated',dim(pm)[1])

vals=as.vector(Y1.real)
n.site=rep(1:(dim(Y1)[1]),nt)
n.time=rep(1:nt,each=(dim(Y1)[1]))
n.type=rep(1,length(as.vector(Y1.real)))
pm1=data.frame(n.site,n.time,n.type,vals) 
pm1$TF=rep('true',dim(pm1)[1])




vals=as.vector(Y2.real)
n.site=rep(1:(dim(Y2)[1]),nt)
n.time=rep(1:nt,each=(dim(Y2)[1]))
n.type=rep(2,length(as.vector(Y2.real)))
pm2=data.frame(n.site,n.time,n.type,vals) 
pm2$TF=rep('true',dim(pm2)[1])

allpm=rbind(pm,pm1,pm2)


# some plots by site and type 
q.sites=sample(1:300,15)
ggplot(allpm %>% filter(n.site%in%q.sites & n.type==1))+geom_point(aes(x=n.time,y=vals,col=TF))+
  theme_bw()+facet_wrap(~n.site)+scale_x_continuous(breaks=seq(1:nt))

q.sites=sample(1:300,15)
ggplot(allpm %>% filter(n.site%in%q.sites & n.type==2))+geom_point(aes(x=n.time,y=vals,col=TF))+
  theme_bw()+facet_wrap(~n.site)+scale_x_continuous(breaks=seq(1:nt))

### all together 

misstime=c()
misssite=c()
misstype=c()
estimated=c()
realval=c()
iteration=c()


for (t in 1:nt)
{
  
  misstime=c(misstime,rep(t,length(miss1[[t]])))
  misstime=c(misstime,rep(t,length(miss2[[t]])))
  misssite=c(misssite,miss1[[t]])
  misssite=c(misssite,miss2[[t]])
  misstype=c(misstype,rep(1,length(miss1[[t]])))
  misstype=c(misstype,rep(2,length(miss2[[t]])))
  
  estimated=c(estimated,keep.Y1.M[j,miss1[[t]],t])
  estimated=c(estimated,keep.Y2.M[j,miss2[[t]],t])
  realval=c(realval,Y1.real[miss1[[t]],t])
  realval=c(realval,Y2.real[miss2[[t]],t])
  
  iteration=c(iteration,rep(j,length(miss1[[t]])+length(miss2[[t]])))
  
  
}


Mdata=data.frame(misstime,misssite,misstype,estimated,realval,iteration)
Mdata$misstype=as.factor(Mdata$misstype)

miss.data=ggplot(Mdata)+geom_point(aes(realval,estimated,col=misstype))+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)+theme_bw()




