library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis)


setwd("/Users/hongjianyang/Research/PA/")
pa2020 <- read.csv('Data/NC_PA_FRM_PM/PA_2020_hourly_PM_NC.csv')
epa2020 <- read.csv('Data/NC_PA_FRM_PM/FRM_2020_hourly_PM_NC.csv')
# Subset data
col <- c('Lon', 'Lat', 'Timestamp', 'PM25_corrected')
pa <- pa2020 %>%
  select(col) %>%
  rename(PM25 = PM25_corrected) %>%
  group_by(Lon, Lat, Timestamp) %>%
  summarise(PM25 = mean(PM25))

pa <- pa[order(pa$Lon, pa$Lat, pa$Timestamp), ]
# Remove NAs
pa <- pa[complete.cases(pa), ]

# Process EPA data
epa2020$Timestamp <- paste(epa2020$Date_GMT, epa2020$Hour_GMT)
epa2020$Timestamp <- paste(epa2020$Timestamp, ':00', sep = '')
col <- c('Lon', 'Lat', 'Timestamp', 'PM25')
epa <- epa2020 %>% 
  select(col)
epa <- epa[order(epa$Lon, epa$Lat, epa$Timestamp), ]
pa$Timestamp <- as.POSIXct(pa$Timestamp)
epa$Timestamp <- as.POSIXct(epa$Timestamp, format = '%Y-%m-%d %H:%M:%OS')

# Prediction locations
# Longitude first, then latitude
s01   <- seq(-90,-65,0.1) # Lon
s02   <- seq(30,38,0.1) # Lat
s0    <- as.matrix(expand.grid(s01,s02))
innc <- map.where(map('state', 'north carolina', fill = TRUE), s0[,1], s0[,2])
s0    <- s0[!is.na(innc),]

X0     <- cbind(1,s0)
X0[,2] <- (X0[,2]-mean(X0[, 2])) / sd(X0[,2])
X0[,3] <- (X0[,3]-mean(X0[,3])) / sd(X0[,3])

lon <- X0[, 2]
lat <- X0[, 3]
lon2 <- lon * lon
lat2 <- lat * lat
lonlat <- lon * lat

X0 <- cbind(1, lat, lon, lat2, lon2, lonlat)

# Get 07/01/2020 - 07/07/2020 data

startTime <- as.POSIXct('2020-07-05 02:00:00')
endTime <- startTime
pa <- subset(pa, Timestamp <= endTime & Timestamp >= startTime)
epa <- subset(epa, Timestamp <= endTime & Timestamp >= startTime)

# Stage 1: Kriging on Purple Air data
S <- data.matrix(pa[, 1:2])
maxd      <- max(dist(S))

# Construct locations
X <- data.frame(cbind(1, S))
X[, 2] <- (X[, 2] - mean(X[, 2])) / sd(X[, 2])
X[, 3] <- (X[, 3] - mean(X[, 3])) / sd(X[, 3])
lon <- X[, 2]
lat <- X[, 3]
lon2 <- lon * lon
lat2 <- lat * lat
lonlat <- lon * lat
X <- cbind(1, lon, lat, lon2, lat2, lonlat)




# Model: Y = \beta_{0} + \beta_{1}*Lon + \beta_{2}*Lat + \beta_{3}*Lon^2 + 
# \beta_{4}*Lat^2 + \beta_{5}*LonLat
Y <- pa$PM25

n.samples <- 25000
burn      <- 5000

starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(Y), "tau.sq"=0.5*var(Y))
tuning    <- list("phi"=1, "sigma.sq"=0.05*var(Y), "tau.sq"=0.05*var(Y))
amcmc     <- list("n.batch"=n.samples/100, "batch.length"=100, "accept.rate"=0.43)
priors    <- list("beta.Norm"=list(rep(0,6), 100*diag(6)),
                  "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                  "sigma.sq.IG"=c(2, 1),
                  "tau.sq.IG"=c(2, 1))

fit_spbayes  <- spLM(Y~X-1, coords=S,
                     starting=starting, tuning=tuning, 
                     priors=priors, cov.model="exponential",
                     amcmc=amcmc, n.samples=n.samples,verbose=FALSE)

pred <- spPredict(fit_spbayes, pred.coords=s0, pred.covars=X0, 
                  start=burn, thin=10, verbose=FALSE)

epa_s0 <- data.matrix(epa[, 1:2])

pred_epa <- spPredict(fit_spbayes, pred.coords=epa_s0, pred.covars=X0, 
                      start=burn, thin=10, verbose=FALSE)

pred <- pred$p.y.predictive.samples
pred_epa <- pred_epa$p.y.predictive.samples

Yhat  <- apply(pred,1,mean)
epa_hat <- apply(pred_epa,1,mean)
pa_std <- apply(pred, 1, sd)

df <- data.frame(long=s0[,1],lat=s0[,2],Y=Yhat)
dfsd <- data.frame(long=s0[,1],lat=s0[,2],Y=pa_std)

ggplot(df, aes(long, lat)) +
  geom_polygon(aes(x=long, y=lat))+
  geom_raster(aes(fill = Y)) + 
  scale_fill_gradientn(colours = viridis(10))+
  coord_fixed() +
  xlab('epa only prediction 07/05 2am')

ggplot(dfsd, aes(long, lat)) +
  geom_polygon(aes(x=long, y=lat))+
  geom_raster(aes(fill = Y)) + 
  scale_fill_gradientn(colours = viridis(10))+
  coord_fixed() +
  xlab('epa only std 07/05 2am')

# Stage 2: Use Predicted values as covariates
# Model: E(Y) = \beta_{0} + \beta_{1} * lon + \beta_{2} * lat + \beta_{3} * lon2
#               + \beta_{4} * lat2 + \beta5 * lonlat + \beta6 * X(s)
Y <- epa$PM25
priors    <- list("beta.Norm"=list(rep(0,7), 100*diag(7)),
                  "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                  "sigma.sq.IG"=c(2, 1),
                  "tau.sq.IG"=c(2, 1))
lon <- epa_s0[, 1]
lat <- epa_s0[, 2]
lon <- (lon - mean(lon)) / sd(lon)
lat <- (lat - mean(lat)) / sd(lat)
lon2 <- lon * lon
lat2 <- lat * lat
lonlat <- lon * lat

X <- cbind(1, lon, lat, lon2, lat2, lonlat, epa_hat)
S <- epa_s0
fit_spbayes_stage2  <- spLM(Y~X-1, coords=S,
                            starting=starting, tuning=tuning, 
                            priors=priors, cov.model="exponential",
                            amcmc=amcmc, n.samples=n.samples,verbose=FALSE)

lon <- s0[, 1]
lat <- s0[, 2]
lon <- (lon - mean(lon)) / sd(lon)
lat <- (lat - mean(lat)) / sd(lat)
lon2 <- lon * lon
lat2 <- lat * lat
lonlat <- lon * lat

X0 <- cbind(1, lon, lat, lon2, lat2, lonlat, Yhat)
pred_stage2 <- spPredict(fit_spbayes_stage2, pred.coords=s0, pred.covars=X0, 
                         start=burn, thin=10, verbose=FALSE)

pred_stage2 <- pred_stage2$p.y.predictive.samples
Yhat <- apply(pred_stage2, 1, mean)

std <- apply(pred_stage2, 1, sd)

df <- data.frame(long=s0[,1],lat=s0[,2],Y=Yhat)
dfsd <- data.frame(long=s0[,1],lat=s0[,2],Y=std)

ggplot(df, aes(long, lat)) +
  geom_polygon(aes(x=long, y=lat))+
  geom_raster(aes(fill = Y)) + 
  scale_fill_gradientn(colours = viridis(10))+
  coord_fixed()

ggplot(dfsd, aes(long, lat)) +
  geom_polygon(aes(x=long, y=lat))+
  geom_raster(aes(fill = Y)) + 
  scale_fill_gradientn(colours = viridis(10))+
  coord_fixed() +
  xlab(startTime)