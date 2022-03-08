## MCMC Gibbs sampling

library(mvtnorm)
library(truncnorm)
setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Spectral/simulate_process.R')


# 1. al
# prior parameters
sigma=5
mu=1
# full conditional distribution
P=t(u2)%*%u2/(tau2)+1/(sigma)
invP=1/P
W=t(Y2)%*%u2/(tau2)-t(v2)%*%u2/(tau2)+mu/(sigma)
invP*W
hist(rtruncnorm(1000, a=0, b=Inf, mean = invP*W, sd = sqrt(invP)))
hist(rnorm(1000,invP*W,sqrt(invP)))
abline(v=al,col='red')


#2. sigmav
# prior parameters
alpha=0.01
beta=0.01
# full conditional distribution
k=length(v2)
p1v=alpha+k/2
invSigmav22=solve(Sigmav22)
p2v=beta+t(v2)%*%invSigmav22%*%v2/2

hist(1/rgamma(1000,p1v,p2v))
abline(v=sigmav,col='red')


#sigmau 
# prior parameters
a=0.1
b=0.1
# full conditional distribution
nk=length(u)
p1u=a+nk/2
invM=solve(M)
p2u=b+t(u)%*%invM%*%u/2

hist(1/rgamma(1000,p1u,p2u))
abline(v=sigmau,col='red')


plot(1/rgamma(1000,p1u,p2u),type='l')
abline(h=sigmau,col='red')


#u1l
Sigmau11=M[1:n[1],1:n[1]]
Sigmau12=M[1:n[1],(n[1]+1):(n[1]+n[2])]
Sigmau22=M[(n[1]+1):(n[1]+n[2]),(n[1]+1):(n[1]+n[2])]

invSigmau22=solve(Sigmau22)
S1=Sigmau11-Sigmau12%*%invSigmau22%*%t(Sigmau12)
Pu1=diag(n[1])/(tau1)+solve(S1)/(sigmau)
Pu1inv=solve(Pu1)

A1=Sigmau12%*%solve(Sigmau22)
Wu1=Y1/(tau1)+t(solve(S1))%*%A1%*%u2/(sigmau)

U1s=rmvnorm(1000,mean=Pu1inv%*%Wu1,sigma=Pu1inv)

par(mar=c(1,1,1,1))
par(mfrow=c(10,10))
for (i in 1:n[1])
{
  hist(U1s[,i])
  abline(v=u1[i],col='red')
}
  

#u2l
S2=Sigmau22-t(Sigmau12)%*%solve(Sigmau11)%*%Sigmau12
A2=t(Sigmau12)%*%solve(Sigmau11)
Pu2=al^2*diag(n[2])/tau2+solve(S2)/(sigmau)
Wu2=al*Y2/(tau2)-al*v2/(tau2)+solve(S2)%*%A2%*%u1/(sigmau)
Pu2inv=solve(Pu2)

U2s=rmvnorm(1000,mean=Pu2inv%*%Wu2,sigma=Pu2inv)

par(mar=c(1,1,1,1))
par(mfrow=c(10,10))
for (i in 1:n[2])
{
  hist(U2s[,i])
  abline(v=u2[i],col='red')
}


#v2l
Pv2=diag(n[2])/(tau2^2)+solve(Sigmav22)/(sigmav^2)
Pv2inv=solve(Pv2)
Wv2=(Y2-al*u2)/(tau2^2)

V2s=rmvnorm(1000,mean=Pv2inv%*%Wv2,sigma=Pv2inv)

par(mar=c(1,1,1,1))
par(mfrow=c(10,10))
for (i in 1:n[2])
{
  hist(V2s[,i])
  abline(v=v2[i],col='red')
}


# tau1
#priors 
alpha1=.1
beta1=.1

p1.1=n[1]/2+alpha1
p2.1=t(Y1-u1)%*%(Y1-u1)/2+beta1

hist(1/rgamma(1000,p1.1,p2.1))
abline(v=tau1,col='red')

#tau2 
#priors
alpha2=.1
beta2=.1


p1.2=n[2]/2+alpha2
p2.2=t(Y2-(al*u2+v2))%*%(Y2-(al*u2+v2))/2+beta2

hist(1/rgamma(1000,p1.2,p2.2))
abline(v=tau2,col='red')


