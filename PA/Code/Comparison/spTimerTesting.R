# Compare different algorithm
# Last Update: 08/23/2022
################################################
###########   System Parameters
################################################
#rm(list  = ls())
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Comparison/")
#source("LMC.R")
dat = load("comparison.RData")
n1 = dim(Y1_train)[1] # training type 1 data
n2 = dim(Y2)[1] # training type 2 data
n3 = dim(Y1_test)[1] # testing data
ts = dim(Y1_train)[2] # Time steps
burn = 1000
df1 = data.frame(matrix(rep(0,16), ncol = 4, nrow = 4))
colnames(df1) <- c('lmc', 'mean', 'twos', 'spTimer')
row.names(df1) = c("ConverageProbability", "Variance", "MSE", 'elapsed')


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
##### Testing with only 1 timestep
n1 = 71
n2 = 70
n3 = 29
ts = 1
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
priors<-spT.priors(model="GP",inv.var.prior=Gamm(2,1),
                   beta.prior=Norm(0,10))

initials<-spT.initials(model="GP", sig2eps=0.1,
                       sig2eta=0.5, beta=NULL)

spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)

# Time data
#time.data <- spT.time(t.series=ts,segment=1)
# Call spTimer
start4 = proc.time()[3]
# Double check here
spt <- spT.Gibbs(formula = y.train ~ lon+lat+lon2+lat2+lonlat+idct, model = "GP",
                 data = df.train, coords = ~lon+lat, cov.fnc="exponential",
                 priors = priors, initials = initials, spatial.decay = spatial.decay,
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
  print(var(obs - truth))
  sumV4 = sumV4 + var(obs - truth)
  mse4 = mse4 + sum((obs - truth)^2) / (n3-1)
}
avgV4 = sumV4 / ts
mse4 = mse4 / ts

df1[1, 4] = cover4; df1[2, 4] = avgV4; df1[3, 4] = mse4; df1[4, 4] = t4

beta_hat = c(2.5503, -1.4758, 3.5358,
             1.2640, 1.2666, 4.5214, 0.8426)
lon = coords1_test[, 1]; lat = lon = coords1_test[, 2]
lon2 = lon * lon; lat2 = lat * lat; lonlat = lon * lat


kappa = function(phi, nu) {
  const = 2
  first = 1/(2^(nu-1) * gamma(nu))
  second = (2 * sqrt(nu) * const * phi)^nu
  third = besselK(2 * sqrt(nu) * const * phi, nu) # K_n(x); n is the second parameter
  return(first * second * third)
}
