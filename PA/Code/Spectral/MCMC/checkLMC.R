# LMC function performs bad. And this script is used to check problems
# First will check kriging codes in the prediction block
rm(list = ls())
library(ggplot2)
library(viridis)
##############################################:
#####           PREDICTIONS            #######:
##############################################:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("simAllTS.R")
# Plot simulated data
real= rbind(Y1.real, Y2.real)
df = data.frame(cbind(coords, real[, 2]))
ggplot(df, aes(X1, X2)) +
  geom_point(aes(colour = X3), size = 0.6) +
  scale_colour_gradient(low="#22FF00", high="#FF0000")
## Bookkeeping
# Prediction coords
s01   <- seq(0,1,0.01) # Lon
s02   <- seq(0,1,0.01) # Lat
sp1 = as.matrix(expand.grid(s01,s02))
nt = 10
n1 = 300
n2 = 350
np1 = dim(sp1)[1]
np2 = 0  # Since we don't want to predict purple air stations
U1p = matrix(0,np1,nt)
Z1p = matrix(0,np1,nt)
Ys1.pred = matrix(0,np1,nt)

all.d = as.matrix(dist(rbind(coords1,coords2,sp1))) 
rangeu=0.9
Mp=exp_corr(all.d, range=rangeu)

Mp00=Mp[1:(n1+n2),1:(n1+n2)]
Mp11=Mp[(n1+n2+1):(n1+n2+np1+np2),(n1+n2+1):(n1+n2+np1+np2)]
Mp10=Mp[(n1+n2+1):(n1+n2+np1+np2),1:(n1+n2)]
Mp01=t(Mp10)

E00=eigen(Mp00)
E00.G=E00$vectors
E00.D=E00$values
Mp00.inv=E00.G%*%diag(1/E00.D)%*%t(E00.G)
AA=Mp10%*%E00.G
B=Mp11-Mp10%*%Mp00.inv%*%Mp01
Uls=rbind(u1,u2)  

for (r in 1:nt) {
  Au=AA%*%Uls[,r] 
  sigmaB=sigmau[r]*B
  Ul.pred=rmvnorm(1,mean=Au,sigma=sigmaB)
  U1p[,r]=Ul.pred[(1:np1)]
}
for(i in 1:nt){Z1p[i,] <- fft_real(U1p[i,],inverse=TRUE)}
Y1.pred <- data.frame(Z1p+rnorm(n=nt,sd=sqrt(tau1)))

ggplot(df, aes(X1, X2)) +
  geom_point(aes(colour = X3), size = 0.6) +
  scale_colour_gradient(low="#22FF00", high="#FF0000") + 
  geom_raster(data = Y1.pred, aes(long, lat), aes(fill = X1)) +
  scale_fill_gradientn(colours = viridis(10))
  
df1 = data.frame(cbind(sp1, Y1.pred[, 1]))
ggplot(df1, aes(s01, s02)) +
  geom_point(aes(colour = V3), size = 0.1) +
  scale_colour_gradient(low="#22FF00", high="#FF0000")


