
#### simulation of the process (just for one time of frequency)

library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)

# set some parameters
n=c(25,15) # number of locations


tau1=1 # error variance1
tau2=1.5  # error variance2

# correlation parameters
sigma1=1 
sigma2=4
rangeu=exp(0.1)
rangev=exp(0.3)

# simulate coordinates
set.seed(22)
coords1 = cbind(runif(n[1],0,1), runif(n[1],0,1))
coords2 = cbind(runif(n[2],0,1), runif(n[2],0,1))
coords=rbind(coords1,coords2)

# Get U and V

# correlation function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}



### 
du1=as.matrix(dist(coords1))
Sigmau11=sigma1*apply(du1,2,exp_corr,range=rangeu)

du2=as.matrix(dist(coords2))
Sigmau22=sigma1*apply(du2,2,exp_corr,range=rangeu)

dv2=as.matrix(dist(coords2))
Sigmav22=sigma1*apply(dv2,2,exp_corr,range=rangev)

u1=mvrnorm(n = 1, rep(0,n[1]), Sigmau11)
#u1 = as.vector(t(chol(Sigmau11)) %*% rnorm(n[1]))
u2=mvrnorm(n = 1, rep(0,n[2]), Sigmau22)
#u2 = as.vector(t(chol(Sigmau22)) %*% rnorm(n[2]))
v2=mvrnorm(n = 1, rep(0,n[2]), Sigmav22)
#v2 = as.vector(t(chol(Sigmav22)) %*% rnorm(n[2]))

# simulate response Y
Y1=u1+rnorm(n[1],tau1)
Y2=u2+v2+rnorm(n[2],tau2)

# set data frame 
type=c(rep('1',n[1]),rep('2',n[2]))
Y=c(Y1,Y2)
U=c(u1,u2)
V=c(rep(NA,n[1]),v2)
s1=c(coords1[,1],coords2[,1])
s2=c(coords1[,2],coords2[,2])

data=data.frame(type,Y,U,V,s1,s2)

# plot data
ggplot(data)+geom_point(aes(x=s1,y=s2,col=Y))+facet_grid(~type)+theme_bw()+ coord_fixed(ratio = 1)


