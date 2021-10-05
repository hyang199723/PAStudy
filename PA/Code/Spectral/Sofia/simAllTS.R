source('ExtraFunctions.R')
#### simulation of the time series process

library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)
library(mvtnorm)

# correlation function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}

###### set some parameters ########

n=c(200,250) # number of locations

nt=10 # total time steps

tau1=1^2 # error variance1
tau2=(1.5)^2  # error variance2
set.seed(99)
al=runif(nt,min=.5,max=1)

# correlation parameters
set.seed(88)
sigmau=rnorm(nt,sd=2)+seq(from=102,to=120,by=2)
set.seed(564)
sigmav=rnorm(nt,sd=2)+seq(from=10,to=28,by=2)
set.seed(213)
rangeu=runif(nt,min=0.5,max=1.5)
lrangeu=log(rangeu)
set.seed(313)
rangev=runif(nt,min=0.5,max=1.5)
lrangev=log(rangev)

## mean
beta.1=0
beta.2=0
###### simulate coordinates #######

set.seed(22)
coords1 = cbind(runif(n[1],0,1), runif(n[1],0,1))
coords2 = cbind(runif(n[2],0,1), runif(n[2],0,1))
coords=rbind(coords1,coords2)

######## Get U and V ##########

du12=as.matrix(dist(coords)) #distance matrix U

u=matrix(NA,ncol=nt,nrow=sum(n))
u1=matrix(NA,ncol=nt,nrow=n[1])
u2=matrix(NA,ncol=nt,nrow=n[2])
v2=matrix(NA,ncol=nt,nrow=n[2])
dv2=as.matrix(dist(coords2)) # distance matrix v

for (t in 1:nt)
{
  #u
  M=exp_corr(du12,range=rangeu[t])
  u[,t]=t(chol(M))%*%rnorm(sum(n),0,sqrt(sigmau[t]))
  u1[,t]=u[(1:n[1]),t]
  u2[,t]=u[(n[1]+1):(sum(n)),t]
  
  #v
  Sigmav22=exp_corr(dv2,range = rangev[t])
  v2[,t]=t(chol(Sigmav22))%*%rnorm(n[2],0,sqrt(sigmav[t]))
  
}



####### simulate response Y ############

#spectral
Z1sp=matrix(NA,ncol=nt,nrow=n[1])
Z2sp=matrix(NA,ncol=nt,nrow=n[2])

for (t in 1:nt)
{
  Z1sp[,t]=u1[,t]#+rnorm(n[1],0,sqrt(nt*tau1))
  Z2sp[,t]=al[t]*u2[,t]+v2[,t]#+rnorm(n[2],0,sqrt(nt*tau2))
  
}

#time domain
Y1=matrix(NA,ncol=nt,nrow=n[1])
Y2=matrix(NA,ncol=nt,nrow=n[2])
for(i in 1:n[1]){Y1[i,] <- fft_real(Z1sp[i,],inverse=TRUE)+beta.1+rnorm(nt,0,sqrt(tau1))}
for(i in 1:n[2]){Y2[i,] <- fft_real(Z2sp[i,],inverse=TRUE)+beta.2+rnorm(nt,0,sqrt(tau2))}

## missing values 

miss1=list()
miss2=list()
set.seed(8923)
for (t in 1:nt)
{
  miss1[[t]]=sample(1:n[1],ceiling(nt/5),replace=FALSE)
  miss2[[t]]=sample(1:n[2],ceiling(nt/5),replace=FALSE)
  Y1[miss1[[t]],t]=NA
  Y2[miss2[[t]],t]=NA
}



###### set data frame and plot #######

time=rep(seq(1,nt),each=sum(n))
type=rep(c(rep('1',n[1]),rep('2',n[2])),nt)
Y=c(as.vector(Y1),as.vector(Y2))
U=c(as.vector(u1),as.vector(u2))
V=c(rep(NA,n[1]*nt),as.vector(v2))
s1=rep(c(coords1[,1],coords2[,1]),nt)
s2=rep(c(coords1[,2],coords2[,2]),nt)
site=rep(seq(1:sum(n)),nt)


data=data.frame(time,type,Y,U,V,s1,s2,site)
s1=coords1
s2=coords2
type=c(rep('1',n[1]),rep('2',n[2]))

# plot some data
#spatial
ggplot(data %>% filter(time==1))+geom_point(aes(x=s1,y=s2,col=Y))+facet_grid(~type)+theme_bw()+ coord_fixed(ratio = 1)
ggplot(data %>% filter(time==2))+geom_point(aes(x=s1,y=s2,col=Y))+facet_grid(~type)+theme_bw()+ coord_fixed(ratio = 1)
ggplot(data %>% filter(time==3))+geom_point(aes(x=s1,y=s2,col=Y))+facet_grid(~type)+theme_bw()+ coord_fixed(ratio = 1)
ggplot(data %>% filter(time==4))+geom_point(aes(x=s1,y=s2,col=Y))+facet_grid(~type)+theme_bw()+ coord_fixed(ratio = 1)

#temporal
ggplot(data %>% filter(site %in%seq(1:20)))+geom_line(aes(x=time,y=Y))+facet_wrap(~site)+theme_bw()




