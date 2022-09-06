rm(list = ls())
# Generate simulation data from different models
## Required data
# coords1: All PA (Type 1) coordinates
# coords1_test: subset of coords1: testing coordinates
# coords1_train: subset of coords1: training coordinates
# coords2: All EPA (Type 2) coordinates
# Y1_sum: All PA (Type 1) data
# Y1_test: Type 1 testing data
# Y1_train: Type 1 training data
# Y2: All Type 2 data

############################################################
#### spTimer AR Models
# Model: Y = O_t + e
# O_t = \rho * (O_{t-1} -Xbeta ) + X*\beta + \eta
# \eta is a spatial process
############################################################
set.seed(123)
a1 = 100
a2 = 70
# number of locations
n=c(a1,a2)
ts = 15
# Randomly select 100 testing locations
vld = rbinom(a1, 1, 0.302)
nt=ts # total time steps
# Error variance
sig2eps=3^2
eps = matrix(rnorm(sum(n) * ts, 0, sqrt(sig2eps)), nrow = sum(n), ncol = ts)
# Spatial
range1=2

# rho
rho = 0.6
# simulation
set.seed(1)
leng = 15
coords1 = cbind(runif(n[1],0,leng), runif(n[1],0,leng))
set.seed(28)
coords2 = cbind(runif(n[2],0,leng), runif(n[2],0,leng))
coords=rbind(coords1,coords2)
# Mean
## mean
lon = c(coords1[, 1], coords2[ ,1])
lat = c(coords1[, 2], coords2[ ,2])
lon2 = lon^2; lat2 = lat^2; lonlat = lon * lat
idct = c(rep(0, a1), rep(1, a2))
X = cbind(rep(1, sum(n)), lon, lat, lon2, lat2, lonlat, idct)
sig2eta = 5^2
set.seed(123)
beta = rnorm(7, 1, 2) / 100
beta[7] = sqrt(sig2eta) / 2
xb = X %*% beta


# Bookkeeping
Y1=matrix(NA,ncol=nt,nrow=n[1])
Y2=matrix(NA,ncol=nt,nrow=n[2])
O_all = matrix(NA,ncol=nt,nrow=sum(n))
# Distance matrix
d = as.matrix(dist(coords))
M = exp(-d / range1)

O_all[, 1] = xb + t(chol(M)) %*% rnorm(sum(n), 0, sqrt(sig2eta))

  for (i in 2:ts) {
    O_previous = O_all[, (i-1)]
    O = rho * (O_previous - xb) + xb + t(chol(M)) %*% rnorm(sum(n), 0, sqrt(sig2eta*(1-rho^2)))
    O_all[, i] = O
  }


# sigma * sqrt(1-rho^2) # make variance the same for all time steps


Y = O_all + eps

Y1 = Y[1:(n[1]), ]; Y2 = Y[(n[1]+1):sum(n), ]
# Train test split
Y1_sum = cbind(Y1, vld)
Y1_train = subset(Y1, vld == 0)
Y1_test = subset(Y1, vld == 1)

c1 = cbind(coords1, vld)
coords1_train = subset(coords1, vld == 0)
coords1_test = subset(coords1, vld == 1)

save(list = c("Y1_sum", "Y1_train", "Y1_test", "Y2", "coords1", "coords1_train", "coords1_test", "coords2"), file = "comparison.RData")

# Plot of simulation data
## Variance by time
v = rep(0, ts)
for (i in 1:ts) {
  v[i] = var(Y[, i])
}
plot(v, main = "AR(1) generated variance")

