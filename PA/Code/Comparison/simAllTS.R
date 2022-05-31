source('ExtraFunctions.R')
#### simulation of the time series process
library(tidyverse)
library(spBayes)
library(ggplot2)
library(mgcv)
library(MASS)
library(mvtnorm)
library(truncnorm)
library(viridis)

# correlation function
exp_corr=function(d,range)
{
  out=exp(-d/range)
  return(out)
}


# Constant correlation across all frequencies.
# Get strong correlation first.

###### set some parameters ########

n=c(200,300) # number of locations

nt=5 # total time steps

tau1=0.5^2 # error variance1
tau2=(0.6)^2  # error variance2
set.seed(99)
#al=runif(nt,min = 5,max=6)
# Change al from uniform sequence to decreasing sequence
al = seq(from = 10, to = 1, length = nt) + rnorm(nt,sd=2)

# correlation parameters
set.seed(88)
sigmau=seq(from=10,to=1,length=nt)+rtruncnorm(nt,a=0,sd=.2)
set.seed(564)
sigmav=seq(from=5,to=.1,length=nt)+rtruncnorm(nt,a=0,sd=.2)
# same range for all freq
rangeu=0.3
lrangeu=log(rangeu)
rangev=0.2
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

M=exp_corr(du12,range=rangeu)
Sigmav22=exp_corr(dv2,range = rangev)
for (t in 1:nt)
{
  #u
  #M=exp_corr(du12,range=rangeu[t])
  u[,t]=t(chol(M))%*%rnorm(sum(n),0,sqrt(sigmau[t]))
  u1[,t]=u[(1:n[1]),t]
  u2[,t]=u[(n[1]+1):(sum(n)),t]
  
  #v
  #Sigmav22=exp_corr(dv2,range = rangev[t])
  v2[,t]=t(chol(Sigmav22))%*%rnorm(n[2],0,sqrt(sigmav[t]))
  
}

####### simulate response Y ############

#spectral
Z1sp=matrix(NA,ncol=nt,nrow=n[1])
Z2sp=matrix(NA,ncol=nt,nrow=n[2])

for (t in 1:nt)
{
  Z1sp[,t]=u1[,t] #+ rnorm(n[1],0,sqrt(nt/2*tau1))
  Z2sp[,t]=al[t]*u2[,t]+v2[,t] #+ rnorm(n[2],0,sqrt(nt/2*tau2))
}

#time domain
Y1=matrix(NA,ncol=nt,nrow=n[1])
Y2=matrix(NA,ncol=nt,nrow=n[2])
for(i in 1:n[1]){Y1[i,] <- fft_real(Z1sp[i,],inverse=TRUE)+beta.1+rnorm(nt,0,sqrt(tau1))}
for(i in 1:n[2]){Y2[i,] <- fft_real(Z2sp[i,],inverse=TRUE)+beta.2+rnorm(nt,0,sqrt(tau2))}

Y1.real=Y1
Y2.real=Y2


## missing values 

# Exclude missing values for testing purpose 
if (F) {
  miss1=list()
  miss2=list()
  set.seed(8923)
  for (t in 1:nt)
  {
    miss1[[t]]=sample(1:n[1],ceiling(n[1]/6),replace=FALSE)
    miss2[[t]]=sample(1:n[2],ceiling(n[2]/6),replace=FALSE)
    Y1[miss1[[t]],t]=NA
    Y2[miss2[[t]],t]=NA
  }
  
}



###### set data frame and plot #######
time = c(rep(1:nt, each = n[1]), rep(1:nt, each = n[2]))
#time=rep(seq(1,nt),each=sum(n))
type=c(rep('1',n[1]*nt),rep('2',n[2]*nt))
Y=c(as.vector(Y1),as.vector(Y2))
Ys=c(as.vector(Z1sp),as.vector(Z2sp))
U=c(as.vector(u1),as.vector(u2))
V=c(rep(NA,n[1]*nt),as.vector(v2))
#s1=rep(c(coords1[,1],coords2[,1]),nt)
#s2=rep(c(coords1[,2],coords2[,2]),nt)
s1 = c(rep(coords1[, 1], nt), rep(coords2[, 1], nt))
s2 = c(rep(coords1[, 2], nt), rep(coords2[, 2], nt))
site=c(rep(seq(1:n[1]),nt), rep((1+n[1]):sum(n),nt))


data=data.frame(time,type,Y,U,V,s1,s2,site, Ys)
s1=coords1
s2=coords2
type=c(rep('1',n[1]),rep('2',n[2]))

# plot some data
#spatial
ggplot(data %>% filter(time==1))+
  geom_point(aes(x=s1,y=s2,col=Y))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1)+
  ggtitle('Y values, time 1')+scale_colour_gradient(low="#22FF00", high="#FF0000")
# Zoom in Type 1 data
ggplot(data %>% filter(time==1 & type == '1'))+
  geom_point(aes(x=s1,y=s2,col=Y))+theme_bw()+coord_fixed(ratio = 1)+ 
  ggtitle('Y values, time 4')+scale_colour_gradient(low="#22FF00", high="#FF0000")
######
ggplot(data %>% filter(time==2))+
  geom_point(aes(x=s1,y=s2,col=Y))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1)+
  ggtitle('Y values, time 2')+scale_colour_gradient(low="#22FF00", high="#FF0000")
# Zoom in Type 1 data
ggplot(data %>% filter(time==2 & type == '1'))+
  geom_point(aes(x=s1,y=s2,col=Y))+theme_bw()+coord_fixed(ratio = 1)+ 
  ggtitle('Y values, time 4')+scale_colour_gradient(low="#22FF00", high="#FF0000")
######
ggplot(data %>% filter(time==3))+
  geom_point(aes(x=s1,y=s2,col=Y))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1)+ 
  ggtitle('Y values, time 3')+scale_colour_gradient(low="#22FF00", high="#FF0000")
# Zoom in Type 1 data
ggplot(data %>% filter(time==3 & type == '1'))+
  geom_point(aes(x=s1,y=s2,col=Y))+theme_bw()+coord_fixed(ratio = 1)+ 
  ggtitle('Y values, time 4')+scale_colour_gradient(low="#22FF00", high="#FF0000")
######
ggplot(data %>% filter(time==4))+
  geom_point(aes(x=s1,y=s2,col=Y))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1)+ 
  ggtitle('Y values, time 4')+scale_colour_gradient(low="#22FF00", high="#FF0000")
# Zoom in
ggplot(data %>% filter(time==4 & type == '1'))+
  geom_point(aes(x=s1,y=s2,col=Y))+theme_bw()+coord_fixed(ratio = 1)+ 
  ggtitle('Y values, time 4')+scale_colour_gradient(low="#22FF00", high="#FF0000")
#temporal
ggplot(data %>% filter(site %in%seq(1:20)))+geom_line(aes(x=time,y=Y))+facet_wrap(~site)+theme_bw()

ggplot(data %>% filter(site %in%seq(1:20)))+geom_line(aes(x=time,y=Ys))+facet_wrap(~site)+theme_bw()

# Plot U data
ggplot(data %>% filter(time==1))+
  geom_point(aes(x=s1,y=s2,col=U))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('U values, time 1') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")

ggplot(data %>% filter(time==2))+
  geom_point(aes(x=s1,y=s2,col=U))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('U values, time 2') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")

# Plot Ys (response value in spectral domain)
ggplot(data %>% filter(time==1))+
  geom_point(aes(x=s1,y=s2,col=Ys))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('Spectral Y values, time 1') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")

# Plot Ys (response value in spectral domain)
ggplot(data %>% filter(time==2))+
  geom_point(aes(x=s1,y=s2,col=Ys))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('Spectral Y values, time 2') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")

ggplot(data %>% filter(time==3))+
  geom_point(aes(x=s1,y=s2,col=Ys))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('Spectral Y values, time 2') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")


ggplot(data %>% filter(time==4))+
  geom_point(aes(x=s1,y=s2,col=Ys))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('Spectral Y values, time 2') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")

ggplot(data %>% filter(time==5))+
  geom_point(aes(x=s1,y=s2,col=Ys))+
  facet_grid(~type)+theme_bw()+coord_fixed(ratio = 1) + ggtitle('Spectral Y values, time 2') +
  scale_colour_gradient(low="#22FF00", high="#FF0000")


test = subset(data, site %in% seq(1,20))
coord = test[, c('s1', 's2')]
distance = as.matrix(dist(coord))
distance = distance[1:20, 1:20]
diag(distance) = 1
min = min(distance)
which(distance == min)



site7 = subset(data, site == 7)




