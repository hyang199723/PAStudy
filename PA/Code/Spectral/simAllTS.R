
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

n=c(600,450) # number of locations

TT=30 # total time steps

tau1=1^2 # error variance1
tau2=(1.5)^2  # error variance2
set.seed(99)
al=runif(TT)

# correlation parameters
set.seed(88)
sigmau=(1/rgamma(TT,1,1))^2
set.seed(564)
sigmav=(1/rgamma(TT,1,1))^2
set.seed(213)
rangeu=runif(TT,min=0.5,max=1.5)
set.seed(313)
rangev=runif(TT,min=0.5,max=1.5)


###### simulate coordinates #######

set.seed(22)
coords1 = cbind(runif(n[1],0,1), runif(n[1],0,1))
coords2 = cbind(runif(n[2],0,1), runif(n[2],0,1))
coords=rbind(coords1,coords2)

######## Get U and V ##########

du12=as.matrix(dist(coords)) #distance matrix U

u=matrix(NA,ncol=TT,nrow=sum(n))
u1=matrix(NA,ncol=TT,nrow=n[1])
u2=matrix(NA,ncol=TT,nrow=n[2])
v2=matrix(NA,ncol=TT,nrow=n[2])
dv2=as.matrix(dist(coords2)) # distance matrix v

for (t in 1:TT)
{
  #u
  M=exp_corr(du12,range=rangeu[t])
  u[,t]=t(chol(M))%*%rnorm(sum(n),0,sqrt(sigmau[t]))
  u1[,t]=u[(1:n[1]),t]
  u2[,t]=u[(n[1]+1):(sum(n)),t]
  
  #v
  Sigmav22=exp_corr(dv2,range = rangev[t])
  v2[,t]=t(chol(Sigmav22))%*%rnorm(n[2],0,sqrt(sigmav))
  
}



####### simulate response Y ############

Y1=matrix(NA,ncol=TT,nrow=n[1])
Y2=matrix(NA,ncol=TT,nrow=n[2])

for (t in 1:TT)
{
  Y1[,t]=u1[,t]+rnorm(n[1],0,sqrt(tau1))
  Y2[,t]=al[t]*u2[,t]+v2[,t]+rnorm(n[2],0,sqrt(tau2))
  
}

## missing values 

miss1=list()
miss2=list()
set.seed(8923)
for (t in 1:TT)
{
  miss1[[t]]=sample(1:n[1],ceiling(TT/5),replace=FALSE)
  miss2[[t]]=sample(1:n[2],ceiling(TT/5),replace=FALSE)
  Y1[miss1[[t]],t]=NA
  Y2[miss2[[t]],t]=NA
}



###### set data frame and plot #######

time=rep(seq(1,TT),each=sum(n))
type=rep(c(rep('1',n[1]),rep('2',n[2])),TT)
Y=c(as.vector(Y1),as.vector(Y2))
U=c(as.vector(u1),as.vector(u2))
V=c(rep(NA,n[1]*TT),as.vector(v2))
s1=rep(c(coords1[,1],coords2[,1]),TT)
s2=rep(c(coords1[,2],coords2[,2]),TT)
site=rep(seq(1:sum(n)),TT)


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




