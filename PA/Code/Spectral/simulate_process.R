
#### simulation of the process (just for one time of frequency)

library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)

                      ###### set some parameters ########

n=c(300,150) # number of locations

tau1=1^2 # error variance1
tau2=(1.5)^2  # error variance2
al=0.8

# correlation parameters
sigmau=3^2
sigmav=2^2
rangeu=0.4
rangev=0.3


                    ###### simulate coordinates #######

set.seed(22)
coords1 = cbind(runif(n[1],0,1), runif(n[1],0,1))
coords2 = cbind(runif(n[2],0,1), runif(n[2],0,1))
coords=rbind(coords1,coords2)

                     ######## Get U and V ##########

# correlation function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}

du12=as.matrix(dist(coords)) #distance matrix U

M=exp_corr(du12,range=rangeu)
u=t(chol(M))%*%rnorm(sum(n),0,sqrt(sigmau))
u1=u[1:n[1]]
u2=u[(n[1]+1):(sum(n))]

dv2=as.matrix(dist(coords2)) # distance matrix v
Sigmav22=exp_corr(dv2,range = rangev)
v2=t(chol(Sigmav22))%*%rnorm(n[2],0,sqrt(sigmav))

                ####### simulate response Y ############

Y1=u1+rnorm(n[1],0,sqrt(tau1))
Y2=al*u2+v2+rnorm(n[2],0,sqrt(tau2))

               ###### set data frame and plot #######

type=c(rep('1',n[1]),rep('2',n[2]))
Y=c(Y1,Y2)
U=c(u1,u2)
V=c(rep(NA,n[1]),v2)
s1=c(coords1[,1],coords2[,1])
s2=c(coords1[,2],coords2[,2])

data=data.frame(type,Y,U,V,s1,s2)

# plot data
ggplot(data)+geom_point(aes(x=s1,y=s2,col=Y))+facet_grid(~type)+theme_bw()+ coord_fixed(ratio = 1)


