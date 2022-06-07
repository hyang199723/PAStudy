source('simAllTS.R')
source('LMC_function.R')

# sp1=NULL
# sp2=NULL
# mean_range=0
# sd_range=1
# mean_var=0
# sd_var=1
# mean_rho=0
# sd_rho=10
# iters=3000
# burn=1000
# thin=1
# update=10
iters=6000
thin=2
out=LMC_fit(Y1,Y2, s1,s2,sp1=NULL,sp2=NULL,
                 mean_range=0, sd_range=1, mean_var=0, sd_var=1, mean_rho=0,
                 sd_rho=10, iters=iters, burn=1000, thin=thin, update=10)

# create df to plot 
RangeU=as.vector(out$rangeU)
RangeV=as.vector(out$rangeV)
Tau1=as.vector(out$tau1)
Tau2=as.vector(out$tau2)
chain=rep(rep(1:thin,each=iters),6)
index=rep(rep(1:iters,thin),6)


df=data.frame(RangeU,RangeV,Tau1,Tau2,chain,index)
df$chain=as.factor(df$chain)

ggplot(df)+geom_line(aes(x=index,RangeU,color=chain))+
  geom_hline(yintercept=rangeu,linetype="dashed", color = "red")+
  theme_bw()


ggplot(df)+geom_line(aes(x=index,RangeV,color=chain))+
  geom_hline(yintercept=rangev,linetype="dashed", color = "red")+
  theme_bw()


ggplot(df)+geom_line(aes(x=index,Tau1,color=chain),alpha=.5)+
  geom_hline(yintercept=tau1,linetype="dashed", color = "red")+
  theme_bw()#+ylim(c(0,.5))


ggplot(df)+geom_line(aes(x=index,Tau2,color=chain),alpha=.5)+
  geom_hline(yintercept=tau2,linetype="dashed", color = "red")+
  theme_bw()#+ylim(c(0,.5))


####### sigmau, sigmav and Al

# SigmaU chains

par(mfrow=c(6,5),mar=c(1,1,1,1))

for (i in 1:nt)
{
  plot(out$sigmaU[i,,1],type='l',col='lightcoral')
  lines(out$sigmaU[i,,2],type='l',col="#66CCCC")
  abline(h=sigmau[i],lty=2,col='black')
  
}

# SigmaV chains
par(mfrow=c(6,5),mar=c(1,1,1,1))

for (i in 1:nt)
{
  plot(out$sigmaV[i,,1],type='l',col='lightcoral')
  lines(out$sigmaV[i,,2],type='l',col="#66CCCC")
  abline(h=sigmav[i],lty=2,col='black')
  
}


# Al chains
par(mfrow=c(6,5),mar=c(1,1,1,1))

for (i in 1:nt)
{
  plot(out$A[i,,1],type='l',col='lightcoral')
  lines(out$A[i,,2],type='l',col="#66CCCC")
  abline(h=al[i],lty=2,col='black')
  
}





####### U and V 


step=seq(1,nt)
location=sample(1:50,25)


par(mfrow=c(5,5),mar=c(1,1,1,1))
for(j in 1:25)
{
  for (i in 1:3)
  {
    plot(out$U1[location[j],step[i],,1],type='l',col='lightcoral')
    lines(out$U1[location[j],step[i],,2],type='l',col="#66CCCC")
    abline(h=u1[location[j],step[i]],lty=2,col='black')
    
  }
}




par(mfrow=c(5,5),mar=c(1,1,1,1))
for(j in 1:25)
{
  for (i in 1:3)
  {
    plot(out$U2[location[j],step[i],,1],type='l',col='lightcoral')
    lines(out$U2[location[j],step[i],,2],type='l',col="#66CCCC")
    abline(h=u2[location[j],step[i]],lty=2,col='black')
    
  }
}

par(mfrow=c(5,5),mar=c(1,1,1,1))
for(j in 1:25)
{
  for (i in 1:3)
  {
    plot(out$V2[location[j],step[i],,1],type='l',col='lightcoral')
    lines(out$V2[location[j],step[i],,2],type='l',col="#66CCCC")
    abline(h=v2[location[j],step[i]],lty=2,col='black')
    
  }
}