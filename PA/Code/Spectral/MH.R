### Metropolis Hasting 

iters=5000
nchains=2

keep.range1=matrix(NA,iters,nchains)
keep.range2=matrix(NA,iters,nchains)


source('simulate_process.R')

# Initial values

Ru = as.vector(u)/sqrt(sigmau)
Rv = as.vector(v2)/sqrt(sigmav)

for (j in 1:nchains)
{
  
  lrange1 = 0
  lrange2 = 0
for(i in 1:iters){
  # range1
  Ms=exp_corr(du12,range=exp(lrange1))
  curll = dmvnorm(Ru,rep(0,sum(n)),Ms,log=TRUE)
  canrange1 = rnorm(1,lrange1,0.5)
  canM = exp_corr(du12,range=exp(canrange1))
  canll = dmvnorm(Ru,rep(0,sum(n)),canM,log=TRUE)

  MH1 <- canll-curll+dnorm(canrange1,log=TRUE)-dnorm(lrange1,log=TRUE)

if (log(runif(1))<MH1)
{
  lrange1=canrange1
}
  keep.range1[i,j]  <-exp(lrange1)

  # range2
  Ss=exp_corr(dv2,range = exp(lrange2))
  curll2 = dmvnorm(Rv,rep(0,n[2]),Ss,log=TRUE)
  canrange2 = rnorm(1,lrange2,0.5)
  canS = exp_corr(dv2,range=exp(canrange2))
  canll2 = dmvnorm(Rv,rep(0,n[2]),canS,log=TRUE)
  
  MH2 <- canll2-curll2+dnorm(canrange2,log=TRUE)-dnorm(lrange2,log=TRUE)
  
  if (log(runif(1))<MH2)
  {
    lrange2=canrange2
  }
  keep.range2[i,j]  <-exp(lrange2)
  
  
  }
}


# plot(keep.range2[,1],type='l')
# lines(keep.range2[,2],type='l',col='blue')
# abline(h=rangev,col='red')

# data frame of the results

prange1=as.vector(keep.range1)
prange2=as.vector(keep.range2)
post=c(prange1,prange2)
chain=c(rep(rep(seq(1:nchains),each=iters),2))
param=c(rep('range1',iters*nchains),rep('range2',iters*nchains))
index=rep(seq(1,iters),nchains*2)


res=data.frame(post=post,chain=chain,parameter=param,index=index,real=rep(c(rangeu,rangev),each=nchains*iters))
res$chain=as.factor(res$chain)

# plot chains

plotchains=ggplot(res)+geom_line(aes(x=index,y=post,col=chain))+facet_grid(~parameter)+theme_bw()+ 
  geom_hline(aes(yintercept = real))

plotchains
#ggsave('plotchains.png',plotchains)
