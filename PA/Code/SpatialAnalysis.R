setwd("/Users/hongjianyang/Research/AirPollution/")
library(ggplot2)
library(maps)
library(viridis)
library(dplyr)
library(geoR)
library(mapproj)
library(usmap)

nc <- read.csv("Data/nc_sensors.csv", header = TRUE)

df <- subset(nc, location_type == 0, select = c('latitude', 'longitude', 
                                                'pm2.5_24hour'))

df <- df %>%
  rename(pm2.5 = pm2.5_60minute)

Y <- df$pm2.5
s <- as.matrix(df[,1:2])
lon    <- s[,1]
lat    <- s[,2]
lon2   <- s[,1]^2
lat2   <- s[,2]^2
lonlat <- s[,1]*s[,2]

ols <- lm(Y~lon+lat+lon2+lat2+lonlat)
res <- Y-ols$fitted.values
##########
# Variogram
#########



L     <- 20
d_max <- 0.3
d     <- seq(0, d_max, length=L)

vg <- variog(coords = s, data = res, uvec = d) 
plot(vg)

################################
# Kriging
################################
epa <- read.csv("Data/EPA.csv", header = TRUE)
df_epa <- epa %>%
  subset(Date == '01/12/2021') %>%
  rename(PM2.5 = Daily.Mean.PM2.5.Concentration, latitude = SITE_LATITUDE,
         longitude = SITE_LONGITUDE) %>%
  select(PM2.5, latitude, longitude)
  
  
# Locations of epa sensors
s0 <- as.matrix(df_epa[2:3])
# Locations of purple air sensor
s <- as.matrix(df[,1:2])
X <- cbind(lon,lat,lon2,lat2,lonlat)

lon    <- s0[,1]
lat    <- s0[,2]
lon2   <- s0[,1]^2
lat2   <- s0[,2]^2
lonlat <- s0[,1]*s0[,2]
X0     <- cbind(lon,lat,lon2,lat2,lonlat)

fit_mle  <- likfit(data=Y,trend= ~X,coords=s,
                   fix.nugget=FALSE,nugget=10,
                   cov.model="matern",
                   ini = c(0.1, 2))


pred <- krige.conv(data=Y,coords=s, # Describe training data
                   locations=s0,    # Describe prediction sites
                   krige=krige.control(trend.d = ~X,  # Covariates at s
                                       trend.l = ~X0, # Covariates at s0
                                       cov.model="exponential",
                                       beta=fit_mle$beta,
                                       cov.pars=fit_mle$cov.pars,
                                       nugget=fit_mle$nugget))


Yhat <- pred$predict         # Kriging predictions
df_epa$pred <- Yhat
df_epa$res <- df_epa$PM2.5 - df_epa$pred
mse <- sum(df_epa$res^2) / dim(df_epa)[1] - 1

ggplot(data = df_epa, aes(longitude, latitude)) +
  borders('state', col = 'blue', fill = 'grey') +
  geom_point(aes(color = abs(res)), size = 1) +
  scale_colour_gradientn(colours = viridis(3)) +
  coord_map(projection = "albers",lat0 = 34, lat1 = 37,xlim=c(-83.5,-76),ylim=c(34,37))


#################################################
# EPA Mean Trend - bias
#################################################

fit_mle  <- likfit(data=Y-7.25,trend= ~X,coords=s,
                   fix.nugget=FALSE,nugget=10,
                   cov.model="matern",
                   ini = c(0.1, 2))


pred <- krige.conv(data=Y-7.25,coords=s, # Describe training data
                   locations=s0,    # Describe prediction sites
                   krige=krige.control(trend.d = ~X,  # Covariates at s
                                       trend.l = ~X0, # Covariates at s0
                                       cov.model="exponential",
                                       beta=fit_mle$beta,
                                       cov.pars=fit_mle$cov.pars,
                                       nugget=fit_mle$nugget))

Yhat <- pred$predict         # Kriging predictions
df_epa$pred <- Yhat
df_epa$res <- df_epa$PM2.5 - df_epa$pred
mse <- sum(df_epa$res^2) / dim(df_epa)[1] - 1

ggplot(data = df_epa, aes(longitude, latitude)) +
  borders('state', col = 'blue', fill = 'grey') +
  geom_point(aes(color = abs(res)), size = 1) +
  scale_colour_gradientn(colours = viridis(3)) +
  coord_map(projection = "albers",lat0 = 34, lat1 = 37,xlim=c(-83.5,-76),ylim=c(34,37))
