rm(list  = ls())
##
###########################
## Attach library spTimer
###########################
library(spTimer)
###########################
## The GP models:
###########################
##
## Model fitting
##
# Read data
data(NYdata)
# MCMC via Gibbs using default choices
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,
                     data=NYdata, model="GP", coords=~Longitude+Latitude,
                     distance.method = 'euclidean',
                     scale.transform="SQRT")
print(post.gp)
# MCMC via Gibbs not using default choices
# Read data
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE)
DataFit<-subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s)
DataValPred<-subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 8)))
# define the time-series
time.data<-spT.time(t.series=60,segment=1)
# hyper-parameters for the prior distributions
priors<-spT.priors(model="GP",inv.var.prior=Gamm(2,1),
                   beta.prior=Norm(0,10^4))
# initial values for the model parameters
initials<-spT.initials(model="GP", sig2eps=0.01,
                       sig2eta=0.5, beta=NULL, phi=0.001)
# input for spatial decay, any one approach from below
#spatial.decay<-spT.decay(distribution="FIXED", value=0.01)
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)
#spatial.decay<-spT.decay(distribution=Unif(0.01,0.02),npoints=5)
# Iterations for the MCMC algorithms
nItr<-5000
# MCMC via Gibbs
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~ cMAXTMP+WDSP+RH,
                     data=DataFit, model="GP", time.data=time.data,
                     coords=~Longitude+Latitude, priors=priors, initials=initials,
                     nItr=nItr, nBurn=0, report=nItr,
                     tol.dist=2, distance.method="geodetic:km",
                     cov.fnc="exponential", scale.transform="SQRT",
                     spatial.decay=spatial.decay)
print(post.gp)
# Summary and plots
summary(post.gp)
summary(post.gp,pack="coda")
plot(post.gp)
plot(post.gp,residuals=TRUE)
coef(post.gp)
confint(post.gp)
terms(post.gp)
formula(post.gp)
model.frame(post.gp)
model.matrix(post.gp)

# Model selection criteria
post.gp$PMCC
##
## Fit and spatially prediction simultaneously
##
# Read data
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE)
DataFit<-subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s)
DataValPred<-subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 8)))
# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
# MCMC via Gibbs will provide output in *.txt format
# from C routine to avoide large data problem in R
set.seed(11)
post.gp.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,
                             data=DataFit, model="GP", coords=coords,
                             newcoords=pred.coords, newdata=DataValPred,
                             scale.transform="SQRT")
print(post.gp.fitpred)
summary(post.gp.fitpred)
coef(post.gp.fitpred)
plot(post.gp.fitpred)
names(post.gp.fitpred)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.gp.fitpred$prediction[,1]))
######################################
## The GP model for sp class data
######################################
# Creating sp class data
library(sp)
data(meuse)
summary(meuse)
coordinates(meuse) <- ~x+y
class(meuse)
out<-spT.Gibbs(formula=zinc~sqrt(dist),data=meuse,
               model="GP", scale.transform="LOG")
summary(out)
# Create a dataset with spacetime class
library(spTimer)
site<-unique(NYdata[,c("Longitude","Latitude")])
library(spacetime)
row.names(site)<-paste("point",1:nrow(site),sep="")
site <- SpatialPoints(site)
ymd<-as.POSIXct(seq(as.Date("2006-07-01"),as.Date("2006-08-31"),by=1))
# introduce class STFDF
newNYdata<-STFDF(sp=site, time=ymd, data=NYdata) # full lattice
class(newNYdata)
out <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,
                 data=newNYdata, model="GP", scale.transform="SQRT")
summary(out)
###########################
## The AR models:
###########################
##
## Model fitting
##
# Read data
data(NYdata)
# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
# MCMC via Gibbs using default choices
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,
                     data=NYdata, model="AR", coords=coords,
                     scale.transform="SQRT")
print(post.ar)
# MCMC via Gibbs not using default choices
# define the time-series
time.data<-spT.time(t.series=62,segment=1)
# hyper-parameters for the prior distributions
priors<-spT.priors(model="AR",inv.var.prior=Gamm(2,1),
                   beta.prior=Norm(0,10^4))
# initial values for the model parameters
initials<-spT.initials(model="AR", sig2eps=0.01,
                       sig2eta=0.5, beta=NULL, phi=0.001)
# Input for spatial decay
#spatial.decay<-spT.decay(distribution="FIXED", value=0.01)
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)
#spatial.decay<-spT.decay(distribution=Unif(0.01,0.02),npoints=5)
# Iterations for the MCMC algorithms
nItr<-5000
# MCMC via Gibbs
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH,
                     data=NYdata, model="AR", time.data=time.data,
                     coords=coords, priors=priors, initials=initials,
                     nItr=nItr, nBurn=0, report=nItr,
                     tol.dist=2, distance.method="geodetic:km",
                     cov.fnc="exponential", scale.transform="SQRT",
                     spatial.decay=spatial.decay)
print(post.ar)
# Summary and plots
summary(post.ar)
plot(post.ar)
# Model selection criteria
post.ar$PMCC
##
## Fit and spatially prediction simultaneously
##
# Read data
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE)
DataFit<-subset(DataFit, with(DataFit, !(Day %in% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s)
DataValPred<-subset(DataValPred, with(DataValPred, !(Day %in% c(30, 31) & Month == 8)))
# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
# MCMC via Gibbs will provide output in *.txt format
# from C routine to avoide large data problem in R
set.seed(11)
post.ar.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,
                             data=DataFit, model="AR", coords=coords,
                             newcoords=pred.coords, newdata=DataValPred,
                             scale.transform="SQRT")
print(post.ar.fitpred)
summary(post.ar.fitpred)
names(post.ar.fitpred)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.ar.fitpred$prediction[,1]))
