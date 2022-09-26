library(truncnorm)
library(mvtnorm)
library(spTimer)
library(tidyverse)

##### Analysis Real data

#source( "./Real_data/process_real_data.R")
load("./Real_data/all_reail_init.RData")
source('./Real_data/LMC_Real.R') # last version here

ts=length(unique(all.dataF$Timestamp))
n1=length(sites.epa)  
n2=length(sites.pa)  

# Create df to save results

K.fold=5
k.tots=round(n1/K.fold)
iters = 400
burn=50

values=rep(NA,k.tots*4*4)
model=rep(c('lmc', 'mean', 'twos', 'spTimer'),k.tots*4)
index=rep(c("ConverageProbability", "Variance", "MSE", 'elapsed'),each=k.tots*4)
k.index=rep(rep(1:k.tots,each=4),4)

df1=data.frame(values,model,index,k.index)

for (k in 1:k.tots)
{
# Define train and test set for each k

epa.test=((k-1)*K.fold+1):(K.fold*k)
epa.train=sites.epa[!(sites.epa%in%epa.test)]

lmcY2 = all.dataF %>% dplyr::filter(Type=='PA') %>% dplyr::select(PM25)
Y2=t(matrix(lmcY2$PM25,nrow=ts,ncol=n2))
coords2=all.dataF %>% dplyr::filter(id%in%sites.pa) %>% dplyr::summarise(lon=unique(Lon),lat=unique(Lat))
coords2=cbind(coords2$lon,coords2$lat)

lmcYtrain =all.dataF %>% dplyr::filter(id%in%epa.train) %>% dplyr::select(PM25)
Y1.train=t(matrix(lmcYtrain$PM25,nrow=ts,ncol=length(epa.train)))
coords.train=all.dataF %>% dplyr::filter(id%in%epa.train) %>% dplyr::summarise(lon=unique(Lon),lat=unique(Lat))
coords.train=cbind(coords.train$lon,coords.train$lat)

lmcYtest = all.dataF %>% dplyr::filter(id%in%epa.test) %>% dplyr::select(PM25)
Y1.test=t(matrix(lmcYtest$PM25,nrow=ts,ncol=length(epa.test)))
coords.test=all.dataF %>% dplyr::filter(id%in%epa.test) %>% dplyr::summarise(lon=unique(Lon),lat=unique(Lat))
coords.test=cbind(coords.test$lon,coords.test$lat)

# Detrend 
for (i in 1:ts) {
  Y1.train[, i] = Y1.train[, i] - mean(Y1.train[, i],na.rm=TRUE)
  Y2[, i] = Y2[, i] - mean(Y2[, i],na.rm=TRUE)
  Y1.test[, i] = Y1.test[, i] - mean(Y1.test[, i],na.rm=TRUE)
}

########################## 1.LMC ############################

# fit and predict
start1=proc.time()[3]
fit = LMC_fit(Y1.train, Y2, coords.train, coords2, sp1 = coords.test,
              iters = iters, thin = 1)

t1 = proc.time()[3] - start1

# Measure error
pred = fit$Y1.p
N = dim(pred)[1]
cover1 = 0
sumV1 = 0

# 95% confidence interval
for (i in 1:N) {
  for (j in 1:ts) {
    if(!is.na(Y1.test[i,j])) # avoid NA values
    {cur = pred[i, j, ]
    true = Y1.test[i,j]
     m = mean(cur); s = sqrt(var(cur))
    low = qnorm(0.025, m, s); high = qnorm(0.975, m, s)
    if (!is.na(true) & true > low & true < high) {cover1 = cover1 + 1}
  }
  }
}
prob1 = cover1 /sum(!is.na(Y1.test))# (N * ts)

pred1 = Matrix::rowMeans(pred, dims = 2)  # Reduce 3-D pred to 2-D

# MSE
mse1 = 0
for (i in 1:ts) {
  if(sum(!is.na(Y1.test[,i]))>0)
   {mse1 = mse1 + sum((pred1[, i] - Y1.test[, i])^2,na.rm = TRUE) / sum(!is.na(Y1.test[,i])) #N
}}
mse1 = mse1 / ts

# Variance
sumV1 = 0
for (i in 1:ts) {
  for (n in 1:N)
  {
   if (!is.na(Y1.test[n, i]))
   {sumV1=sumV1+var(pred[n,i,]-Y1.test[n, i])
  }
   }
}
avgV1 = sumV1 / sum(!is.na(Y1.test))

# save values
df1=df1 %>%
  mutate(values=replace(values, model=='lmc' & index=='MSE' & k.index==k, mse1)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='lmc' & index=='ConverageProbability' & k.index==k, prob1)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='lmc' & index=='Variance' & k.index==k, avgV1)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='lmc' & index=='elapsed' & k.index==k, t1)) %>%
  as.data.frame()

rm(fit, pred1,lmcYtrain,lmcYtest,lmcY2,pred,mse1,prob1,sumV1,avgV1,t1)
gc()
print('lmc')
###################### 2.Mean-adjusting algorithm #######################

cover2 = 0
sumV2 = 0
mse2 = 0
start2 = proc.time()[3]

for (i in 1:ts) {
  # Observation for one time stamp
  y1.train.obs = as.vector(Y1.train[, i])
  y1.test.obs = as.vector(Y1.test[, i])
  y2.obs = as.vector(Y2[, i])
  y.train = c(y1.train.obs, y2.obs)
  y.test = y1.test.obs
  
  ## Covariates
  # Training
  idct = c(rep(0, n1-K.fold), rep(1, n2))
  c.train = rbind(coords.train, coords2)
  lon = c.train[, 1]; lat = c.train[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.train = cbind(rep(1,n1-K.fold+n2), lon, lat, lon2, lat2, lonlat, idct)
  # Testing
  c.test = coords.test
  idct = rep(0, K.fold)
  lon = c.test[, 1]; lat = c.test[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.test = cbind(rep(1, K.fold), lon, lat, lon2, lat2, lonlat, idct)
  
  #### Chech NA values 
  which.NA=is.na(y.train)
  y.train2=y.train[!which.NA]
  X.train2=X.train[!which.NA,]
  
  which.NA.test=is.na(y.test)
  y.test2=y.test[!which.NA.test]
  
  if (length(y.test2)>0) # avoiding all NA values 
  {
  # Fitting mode
  maxd = max(dist(X.train2))
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(y.train2), "tau.sq"=0.5*var(y.train2))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(y.train2), "tau.sq"=0.05*var(y.train2))
  amcmc     <- list("n.batch"=250, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,7), 100*diag(7)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  
  meanadj  <- spLM(y.train2~X.train2-1, coords=c.train[!which.NA,],
                   starting=starting, tuning=tuning, 
                   priors=priors, cov.model="exponential",
                   amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
  # Predict

  sppred <- spPredict(meanadj, pred.coords=c.test, pred.covars=X.test, 
                      start=burn, thin=10, verbose=FALSE)
  
  sppred <- sppred$p.y.predictive.samples
  Yhat  <- apply(sppred,1,mean)
  Ysd   <- apply(sppred,1,sd)
  # 95% coverage probability and avg varianve
  low = qnorm(0.025, Yhat, Ysd); high = qnorm(0.975, Yhat, Ysd)
  cover2 = cover2 + sum(y.test2 > low[!which.NA.test] & y.test2 < high[!which.NA.test])
  
  # MSE
  mse2 = mse2 + (sum((Yhat[!which.NA.test] - y.test2)^2) / length(y.test2))
  
  for (l in 1:length(y.test2))
  {
  sumV2 = sumV2 + var((sppred[,!which.NA.test]-y.test2)[l,])
  }
  rm(sppred, meanadj)
  gc()
  }
  rm(y1.train.obs, y1.test.obs, y2.obs, y.train, y.test )
  gc()
  
}

t2 = proc.time()[3] - start2
prob2 = cover2 / sum(!is.na(Y1.test))#/(K.fold * ts)
avgV2 = sumV2 / sum(!is.na(Y1.test))#/ts
mse2 = mse2 / ts

# save values
df1=df1 %>%
  mutate(values=replace(values, model=='mean' & index=='MSE' & k.index==k, mse2)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='mean' & index=='ConverageProbability' & k.index==k, prob2)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='mean' & index=='Variance' & k.index==k, avgV2)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='mean' & index=='elapsed' & k.index==k, t2)) %>%
  as.data.frame()

rm(mse2,prob2,avgV2,t2)
print('mean')

######################## 3.Two-stage algorithm #######################

cover3 = 0
sumV3 = 0
mse3 = 0
start3 = proc.time()[3]

for (i in 1:ts) {
  ## Stage 1
  y2.obs = as.vector(Y2[, i])
  y1 = y2.obs
  c1 = coords2
  cp = rbind(coords.train, coords.test)
  ## Covariates
  lon = c1[, 1]; lat = c1[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X1 = cbind(rep(1,n2), lon, lat, lon2, lat2, lonlat)
  # prediction
  lon = cp[, 1]; lat = cp[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  Xp = cbind(rep(1,n1), lon, lat, lon2, lat2, lonlat)
  
  #### Chech NA values 
  which.NA=is.na(y1)
  y1F=y1[!which.NA]
  X1F=X1[!which.NA,]
  
  # Fitting model
  maxd = max(dist(c1[!which.NA,]))
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(y1F), "tau.sq"=0.5*var(y1F))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(y1F), "tau.sq"=0.05*var(y1F))
  amcmc     <- list("n.batch"=250, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,6), 100*diag(6)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  twoS1  <- spLM(y1F~X1F-1, coords=c1[!which.NA,],
                 starting=starting, tuning=tuning, 
                 priors=priors, cov.model="exponential",
                 amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
  # Predict
  s1 <- spPredict(twoS1, pred.coords=cp, pred.covars=Xp, 
                  start=burn, thin=10, verbose=FALSE)
  
  s1 <- s1$p.y.predictive.samples
  Yhat  <- apply(s1,1,mean)
  Yhat.train = Yhat[epa.train]
  Yhat.test = Yhat[epa.test]
  
  # Stage 2
  y1.train.obs = as.vector(Y1.train[, i])
  y1.test.obs = as.vector(Y1.test[, i])
  c.train = coords.train
  # Covariates
  lon = c.train[, 1]; lat = c.train[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.train = cbind(rep(1,n1-K.fold), lon, lat, lon2, lat2, lonlat, Yhat.train)
  
  c.test = coords.test
  lon = c.test[, 1]; lat = c.test[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.test = cbind(rep(1,K.fold), lon, lat, lon2, lat2, lonlat, Yhat.test)
  
  #### Chech NA values 
  which.NA=is.na(y1.train.obs)
  y1.train.obsF=y1.train.obs[!which.NA]
  X.trainF=X.train[!which.NA,]
  
  which.NA.test=is.na(y1.test.obs)
  y1.test.obsF=y1.test.obs[!which.NA.test]
  
  if (length(y1.test.obsF)>0) # avoiding all NA values 
  {
  # Fitting model
  maxd = max(dist(c.train))
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(y1.train.obsF), "tau.sq"=0.5*var(y1.train.obsF))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(y1.train.obsF), "tau.sq"=0.05*var(y1.train.obsF))
  amcmc     <- list("n.batch"=250, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,7), 100*diag(7)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  twoS2  <- spLM(y1.train.obsF~X.trainF-1, coords=c.train[!which.NA,],
                 starting=starting, tuning=tuning, 
                 priors=priors, cov.model="exponential",
                 amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
  # Predict
  s2 <- spPredict(twoS2, pred.coords=c.test, pred.covars=X.test, 
                  start=burn, thin=10, verbose=FALSE)
  
  s2 <- s2$p.y.predictive.samples
  Yhat  <- apply(s2,1,mean)
  Ysd <- apply(s2,1,sd)
  
  # 95% coverage probability and avg varianve
  low = qnorm(0.025, Yhat, Ysd); high = qnorm(0.975, Yhat, Ysd)
  cover3 = cover3 + sum(y1.test.obsF > low[!which.NA.test] & y1.test.obsF < high[!which.NA.test])
  
  # MSE
  mse3 = mse3 + (sum((Yhat[!which.NA.test] - y1.test.obsF)^2) / length(y1.test.obsF))
  
  # Variance
  for (l in 1:length(y1.test.obsF))
  {
    sumV3 = sumV3 + var((s2[,!which.NA.test]-y1.test.obsF)[l,])
  }
  rm(twoS1, s1,twoS2,s2)
  gc()
  }
  #print(i/ts)
}

t3 = proc.time()[3] - start3
prob3 = cover3 / sum(!is.na(Y1.test))#(n3 * ts)
avgV3 = sumV3 / sum(!is.na(Y1.test))#/ts
mse3 = mse3 / ts

# save values
df1=df1 %>%
  mutate(values=replace(values, model=='twos' & index=='MSE' & k.index==k, mse3)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='twos' & index=='ConverageProbability' & k.index==k, prob3)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='twos' & index=='Variance' & k.index==k, avgV3)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='twos' & index=='elapsed' & k.index==k, t3)) %>%
  as.data.frame()

rm(mse3,prob3,avgV3,t3)
print('twos')
######################## 4.spTimer, GP model #######################


# The data should be ordered first by the time
y1.train = as.vector((Y1.train))
y2.train = as.vector((Y2))
y.test = as.vector((Y1.test))
y.train = c(y1.train, y2.train)

## Covariates
# Training
idct = c(rep(0, (n1-K.fold)*ts), rep(1, n2*ts))
c.train = rbind(coords.train, coords2)
lon = c.train[, 1]; lat = c.train[, 2]

# Make duplication for time
lon = rep(lon, each = ts); lat = rep(lat, each = ts)
lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
X.train = cbind(rep(1,(n1-K.fold+n2)*ts), lon, lat, lon2, lat2, lonlat, idct)
df.train = data.frame(cbind(y.train, X.train))
train.coords = rbind(coords.train, coords2)
# Testing
c.test = coords.test
idct = rep(0, K.fold*ts)
lon = c.test[, 1]; lat = c.test[, 2]
lon = rep(lon, each = ts); lat = rep(lat, each = ts)
lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
X.test = data.frame(cbind(rep(1, K.fold*ts), lon, lat, lon2, lat2, lonlat, idct))
df.test = cbind(y.test, X.test)

# Give a more informative prior
priors<-spT.priors(model="AR",inv.var.prior=Gamm(2,1),
                   beta.prior=Norm(0,10))
initials<-spT.initials(model="AR", sig2eps=0.1,
                       sig2eta=0.5, beta=NULL)
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)

start4 = proc.time()[3]
spt <- spT.Gibbs(formula = y.train ~ lon+lat+lon2+lat2+lonlat+idct, model = "AR",
                 data = df.train, coords = ~lon+lat, cov.fnc="exponential",
                 priors = priors, initials = initials, nItr = iters, nBurn = burn,
                 spatial.decay = spatial.decay,
                 distance.method = 'euclidean',tol.dist = 0.0005)


pp=predict(spt,newdata=X.test,newcoords = coords.test)

t4 = proc.time()[3] - start4
preds=pp$pred.samples

mean=apply(preds, 1, mean)
sd=apply(preds, 1, sd)

# Book keeping
# 95% coverage probability
nonas=!is.na(y.test)
low = qnorm(0.025, mean, sd); high = qnorm(0.975, mean, sd)
cover4 = sum(y.test[nonas] > low[nonas] & y.test[nonas] < high[nonas]) / (K.fold*ts)

## MSE

mse4=0
N = K.fold * ts
for (i in 1:ts) {
  q=seq((i-1)*K.fold+1,i*K.fold)
  obs = mean[q]
  truth = y.test[q]
  # obs = mean[seq(i, N, ts)]
  # truth = y.test[seq(i, N, ts)]
  nonas=!is.na(truth)
  mse4 = mse4 + sum((obs[nonas] - truth[nonas])^2) / (length(truth[nonas]))
  print(i)
  print(mse4)
  
}
mse4 = mse4 / ts

## Variance
sumV4 = 0
for (i in 1:length(y.test))
  {
    if (!is.na(y.test[i]))
    {sumV4=sumV4+var(preds[i,]-y.test[i])}
  }

avgV4 = sumV4 / sum(!is.na(y.test))

# save values
df1=df1 %>%
  mutate(values=replace(values, model=='spTimer' & index=='MSE' & k.index==k, mse4)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='spTimer' & index=='ConverageProbability' & k.index==k, cover4)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='spTimer' & index=='Variance' & k.index==k, avgV4)) %>%
  as.data.frame()

df1=df1 %>%
  mutate(values=replace(values, model=='spTimer' & index=='elapsed' & k.index==k, t4)) %>%
  as.data.frame()


rm(pp,mean,sd,preds,spt,t4,avgV4,cover4,mse4)
gc()

print('sptimer')
print(k)

}


#### Plots

ggplot(df1%>% filter(model=='lmc' & index=='Variance'| model=='mean' & index=='Variance'|
                       model=='twos' & index=='Variance',k.index<=6))+
  geom_boxplot(aes(y=values,x=index))+facet_wrap(~model)+
  theme_bw()+ ggtitle('Results Variance' )

ggplot(df1%>% filter(model=='lmc' & index=='MSE'| model=='mean' & index=='MSE'
                     | model=='twos' & index=='MSE',k.index<=6))+
  geom_boxplot(aes(y=values,x=index))+facet_wrap(~model)+
  theme_bw()+ ggtitle('Results MSE' )

ggplot(df1%>% filter(model=='lmc' & index=='elapsed'| model=='mean' & index=='elapsed',k.index<=6))+
  geom_boxplot(aes(y=values,x=index))+facet_wrap(~model)+
  theme_bw()+ ggtitle('Results elapsed' )

ggplot(df1)+ geom_boxplot(aes(y=values,x=index,col=model))+facet_wrap(~index,scales="free")+
  theme_bw()+ ggtitle('All results' )
