# Compare different algorithm
# Last Update: 08/23/2022
################################################
###########   System Parameters
################################################
rm(list  = ls())
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Comparison/")
source("LMC.R")
dat = load("comparison.RData")
n1 = dim(Y1_train)[1] # training type 1 data
n2 = dim(Y2)[1] # training type 2 data
n3 = dim(Y1_test)[1] # testing data
ts = dim(Y1_train)[2] # Time steps
df1 = data.frame(matrix(rep(0,16), ncol = 4, nrow = 4))
colnames(df1) <- c('lmc', 'mean', 'twos', 'spTimer')
row.names(df1) = c("ConverageProbability", "Variance", "MSE", 'elapsed')

################################################
###########   Simulation Parameters
################################################
# Number of MCMC iterations: 7000
iters = 7000
# Burn -in: 2000
burn = 2000
################################################
###########   Start of Comparison Codes
################################################

################################################
###########   LMC algorithm
################################################
notLMC = T
lmcYtrain = Y1_train
lmcY2 = Y2
lmcYtest = Y1_test
# If the simulation data is not from LMC model, detrend

for (i in 1:ts) {
  lmcYtrain[, i] = lmcYtrain[, i] - mean(lmcYtrain[, i])
  lmcY2[, i] = lmcY2[, i] - mean(lmcY2[, i])
  lmcYtest[, i] = lmcYtest[, i] - mean(lmcYtest[, i])
}

start1=proc.time()[3]
fit = LMC_fit(lmcYtrain, lmcY2, coords1_train, coords2, sp1 = coords1_test, 
              iters = iters, thin = 1, burn = burn)
t1 = proc.time()[3] - start1
pred = fit$Y1.p
N = dim(pred)[1]
# 95% prediction coverage and average prediction variance
cover1 = 0
sumV1 = 0

for (i in 1:N) {
  for (j in 1:ts) {
    cur = pred[i, j, ]
    true = lmcYtest[i,j]
    # 95% confidence interval
    m = mean(cur); s = sqrt(var(cur))
    low = qnorm(0.025, m, s); high = qnorm(0.975, m, s)
    if (true > low & true < high) {cover1 = cover1 + 1}
  }
}
prob1 = cover1 / (N * ts)

# Var and mse
# pred = rowMeans(pred, dims = 2)  # Reduce 3-D pred to 2-D
sumV1 = 0; mse1 = 0
for (i in 1:ts) {
  for (j in 1:n3) {
    sumV1 = sumV1 + var(pred[j, i, ]) #  - lmcYtest[j, i]
    mse1 = mse1 + sum((pred[j, i, ] - lmcYtest[j, i])^2) / (iters-burn+1)
  }
}
avgV1 = sumV1 / (ts*n3)
mse1 = mse1 / (ts*n3)

df1[1, 'lmc'] = prob1; df1[2, 'lmc'] = avgV1; df1[3, 'lmc'] = mse1; df1[4, 'lmc'] = t1

# Convergence Diagnostic
range1_fit = fit$rangeU; range2_fit = fit$rangeV; sU = fit$sigmaU; sV = fit$sigmaV
al = fit$A; u1 = fit$U1; u2 = fit$U2; v = fit$V2
sigmaU = sU[,,1]; sigmaV = sV[,,1]; Al = al[,,1]
mean_al = rowMeans(Al)
plot(range1_fit, type = 'l')
plot(range2_fit, type = 'l')
plot(sigmaU[1, ], type = 'l')
plot(sigmaU[7, ], type = 'l')
plot(sigmaV[1, ], type = 'l')
plot(sigmaV[7, ], type = 'l')
plot(Al[1, ], type = 'l')
plot(Al[7, ], type = 'l')
plot(mean_al, type = 'l', main = 'Estimated Al by spectrum')



################################################
###########   Mean-adjusting algorithm
################################################
# 0 means Type 1 station, 1 means Type 2
# Bookkeeping
cover2 = 0
sumV2 = 0
mse2 = 0
start2 = proc.time()[3]
for (i in 1:ts) {
  # Observation for one time stamp
  y1.train.obs = as.vector(Y1_train[, i])
  y1.test.obs = as.vector(Y1_test[, i])
  y2.obs = as.vector(Y2[, i])
  y.train = c(y1.train.obs, y2.obs)
  y.test = y1.test.obs
  ## Covariates
  # Training
  idct = c(rep(0, n1), rep(1, n2))
  c.train = rbind(coords1_train, coords2)
  lon = c.train[, 1]; lat = c.train[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.train = cbind(rep(1,n1+n2), lon, lat, lon2, lat2, lonlat, idct)
  # Testing
  c.test = coords1_test
  idct = rep(0, n3)
  lon = c.test[, 1]; lat = c.test[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.test = cbind(rep(1, n3), lon, lat, lon2, lat2, lonlat, idct)
  # Fitting model
  maxd = max(dist(X.train))
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(y.train), "tau.sq"=0.5*var(y.train))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(y.train), "tau.sq"=0.05*var(y.train))
  # amcmc     <- list("n.batch"=250, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,7), 100*diag(7)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  meanadj  <- spLM(y.train~X.train-1, coords=c.train,
                       starting=starting, tuning=tuning, 
                       priors=priors, cov.model="exponential",
                       n.samples=iters,verbose=FALSE)
  # Predict
  sppred <- spPredict(meanadj, pred.coords=c.test, pred.covars=X.test, 
                    start=burn, thin=1, verbose=FALSE)
  
  sppred <- sppred$p.y.predictive.samples
  Yhat  <- apply(sppred,1,mean)
  YV <- apply(sppred,1,var)
  Ysd   <- apply(sppred,1,sd)
  # 95% coverage probability and avg varianve
  low = qnorm(0.025, Yhat, Ysd); high = qnorm(0.975, Yhat, Ysd)
  cover2 = cover2 + sum(y.test > low & y.test < high)
  
  sumV2 = sumV2 + sum(YV) / n3
  # MSE
  foo = sweep(sppred, 1, y.test, "-")
  foo2 = foo^2
  mse2 = mse2 + sum(rowSums(foo2, 1) / (iters - burn + 1)) / n3
}
t2 = proc.time()[3] - start2
prob2 = cover2 / (n3 * ts)
avgV2 = sumV2 / ts
mse2 = mse2 / ts
df1[1, 'mean'] = prob2; df1[2, 'mean'] = avgV2; df1[3, 'mean'] = mse2; df1[4, 'mean'] = t2
################################################
###########   Two-stage algorithm
################################################
cover3 = 0
sumV3 = 0
mse3 = 0
start3 = proc.time()[3]
for (i in 1:ts) {
  # Stage 1
  y2.obs = as.vector(Y2[, i])
  y1 = y2.obs
  c1 = coords2
  cp = rbind(coords1_train, coords1_test)
  ## Covariates
  lon = c1[, 1]; lat = c1[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X1 = cbind(rep(1,n2), lon, lat, lon2, lat2, lonlat)
  # prediction
  lon = cp[, 1]; lat = cp[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  Xp = cbind(rep(1,n1+n3), lon, lat, lon2, lat2, lonlat)
  # Fitting model
  maxd = max(dist(c1))
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(y1), "tau.sq"=0.5*var(y1))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(y1), "tau.sq"=0.05*var(y1))
  # amcmc     <- list("n.batch"=250, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,6), 100*diag(6)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  twoS1  <- spLM(y1~X1-1, coords=c1,
                       starting=starting, tuning=tuning, 
                       priors=priors, cov.model="exponential",
                       n.samples=iters,verbose=FALSE)
  # Predict
  s1 <- spPredict(twoS1, pred.coords=cp, pred.covars=Xp, 
                      start=burn, thin=1, verbose=FALSE)
  
  s1 <- s1$p.y.predictive.samples
  Yhat  <- apply(s1,1,mean)
  Yhat.train = Yhat[1:n1]
  Yhat.test = Yhat[(n1+1):(n1+n3)]
  
  # Stage 2
  y1.train.obs = as.vector(Y1_train[, i])
  y1.test.obs = as.vector(Y1_test[, i])
  c.train = coords1_train
  # Covariates
  lon = c.train[, 1]; lat = c.train[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.train = cbind(rep(1,n1), lon, lat, lon2, lat2, lonlat, Yhat.train)
  
  c.test = coords1_test
  lon = c.test[, 1]; lat = c.test[, 2]
  lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
  X.test = cbind(rep(1,n3), lon, lat, lon2, lat2, lonlat, Yhat.test)
  
  # Fitting model
  maxd = max(dist(c.train))
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(y1.train.obs), "tau.sq"=0.5*var(y1.train.obs))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(y1.train.obs), "tau.sq"=0.05*var(y1.train.obs))
  # amcmc     <- list("n.batch"=250, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,7), 100*diag(7)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  twoS2  <- spLM(y1.train.obs~X.train-1, coords=c.train,
                       starting=starting, tuning=tuning, 
                       priors=priors, cov.model="exponential",
                  n.samples=iters,verbose=FALSE)
  # Predict
  s2 <- spPredict(twoS2, pred.coords=c.test, pred.covars=X.test, 
                  start=burn, thin=10, verbose=FALSE)
  
  s2 <- s2$p.y.predictive.samples
  Yhat  <- apply(s2,1,mean)
  Ysd <- apply(s2,1,sd)
  YV <- apply(s2, 1, var)
  # 95% coverage probability and avg varianve
  low = qnorm(0.025, Yhat, Ysd); high = qnorm(0.975, Yhat, Ysd)
  cover3 = cover3 + sum(y1.test.obs > low & y1.test.obs < high)
  sumV3 = sumV3 + sum(YV) / n3
  # MSE
  foo = sweep(s2, 1, Yhat, "-")
  foo2 = foo^2
  mse3 = mse3 + sum(rowSums(foo2, 1) / (iters - burn + 1)) / n3
}
t3 = proc.time()[3] - start3
prob3 = cover3 / (n3 * ts)
avgV3 = sumV3 / ts
mse3 = mse3 / ts
df1[1, 'twos'] = prob3; df1[2, 'twos'] = avgV3; df1[3, 'twos'] = mse3; df1[4, 'twos'] = t3

################################################
###########   spTimer, GP model
################################################
library(spTimer)
#The data should be ordered first by the time and then by the sites specified by the coords below. 
# One can also supply coordi- nates through this argument, 
# where coordinate names should be "Latitude" and "Longitude".

# n1: training type 1 data
# n2: training type 2 data
# n3: testing data

# The data should be ordered first by the time
y1.train = as.vector(t(Y1_train))
y2.train = as.vector(t(Y2))
y.test = as.vector(t(Y1_test))
#test = as.vector(Y1_test)
y.train = c(y1.train, y2.train)
## Covariates
# Training
idct = c(rep(0, n1*ts), rep(1, n2*ts))
c.train = rbind(coords1_train, coords2)
lon = c.train[, 1]; lat = c.train[, 2]

# Make duplication for time
lon = rep(lon, each = ts); lat = rep(lat, each = ts)
lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
X.train = cbind(rep(1,(n1+n2)*ts), lon, lat, lon2, lat2, lonlat, idct)
df.train = data.frame(cbind(y.train, X.train))
train.coords = rbind(coords1_train, coords2)
# Testing
c.test = coords1_test
idct = rep(0, n3*ts)
lon = c.test[, 1]; lat = c.test[, 2]
lon = rep(lon, each = ts); lat = rep(lat, each = ts)
lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat
X.test = data.frame(cbind(rep(1, n3*ts), lon, lat, lon2, lat2, lonlat, idct))
df.test = cbind(y.test, X.test)


# Give a more informative prior
priors<-spT.priors(model="AR",inv.var.prior=Gamm(2,1),
                   beta.prior=Norm(0,10))

initials<-spT.initials(model="AR", sig2eps=0.1,
                       sig2eta=0.5, beta=NULL)

spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)

# Time data
#time.data <- spT.time(t.series=ts,segment=1)
# Call spTimer
start4 = proc.time()[3]
# Double check here
spt <- spT.Gibbs(formula = y.train ~ lon+lat+lon2+lat2+lonlat+idct, model = "AR",
                 data = df.train, coords = ~lon+lat, cov.fnc="exponential",
                 priors = priors, initials = initials, nItr = iters, nBurn = burn,
                 spatial.decay = spatial.decay,
                 distance.method = 'euclidean',
                 newcoords=coords1_test, newdata=X.test)
t4 = proc.time()[3] - start4
mean = spt$prediction$Mean
sd = spt$prediction$SD # Check SD here
# Book keeping
cover4 = 0
sumV4 = 0
mse4 = 0
# 95% coverage probability
low = qnorm(0.025, mean, sd); high = qnorm(0.975, mean, sd)
cover4 = sum(y.test > low & y.test < high) / (n3*ts)
# Average variance and MSE
N = n3 * ts
for (i in 1:ts) {
  obs = mean[seq(i, N, ts)]
  truth = y.test[seq(i, N, ts)]

  sumV4 = sumV4 + var(obs - truth)
  mse4 = mse4 + sum((obs - truth)^2) / (n3-1)
}
avgV4 = sumV4 / ts
mse4 = mse4 / ts

df1[1, 4] = cover4; df1[2, 4] = avgV4; df1[3, 4] = mse4; df1[4, 4] = t4

#X = as.matrix(X.test, nrow = 435)
#beta.out = matrix(c(-0.0825, 0.6167, 4.1731, 1.1514, 1.2797, 4.4682, 1.9471), nrow = 7)
#t1 = X %*% beta.out

#a = predict.spT(spt, X.test, coords1_test)

library(xtable)
xtable(df1)
