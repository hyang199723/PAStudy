library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis)

setwd("/Users/hongjianyang/Research/PA/")
source('Code/Tools/myFunc.R')
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



# Check for high variance pa monitors
# Pivot epa wider
a <- data.frame(pivot_wider(pa, id_cols = Timestamp, names_from = c(Lon, Lat), values_from = PM25))
a <- a[order(a$Timestamp), ]
# Get some statistics on a
a.mean <- colMeans(a[, 2:dim(a)[2]], na.rm = TRUE)
a.NA <- colSums(!is.na.data.frame(a))
sum(a.NA == 1) # 18 columns only have 1 obs
# Remove all columns with less than 10 obs
colIndex <- names(a.NA[a.NA <= 10])
a <- a[,!(names(a) %in% colIndex)]

sdNA <- function(col) {
  return(sd(col, na.rm = TRUE))
}

a.sd <- as.vector(apply(a[, 2:dim(a)[2]], 2, sdNA))
# Standard deviation summary:
# first quantile: 3.27
# third quantile: 4.576
# IQR = 1.3
# Any standard deviation above 7 would be considered an outlier
outlier.index <- which(as.numeric(a.sd > 7) %in% c(1)) + 1
a <- a[-outlier.index]

# After removing stations with only 1 observations and outliers, 116 stations left
a.back <- a %>%
  pivot_longer(!Timestamp, names_to = 'loc', values_to = 'PM25')

a.back <- a.back[complete.cases(a.back), ]
temp <- str_remove(a.back$loc, 'X.')
temp <- as.numeric(unlist(strsplit(temp, split = '_')))
a.back$Lon <- temp[c(TRUE, FALSE)]
a.back$Lat <- temp[c(FALSE, TRUE)]
a <- a.back %>%
  select(-c('loc')) %>%
  select(c('Lon', 'Lat', 'Timestamp', 'PM25')) %>%
  mutate(Lon = -Lon)

############################################################
############################################################
# Look at the variance again

pa <- a
epa2020 <- read.csv('Data/NC_PA_FRM_PM/FRM_2020_hourly_PM_NC.csv')
# Subset data
col <- c('Lon', 'Lat', 'Timestamp', 'PM25_corrected')

# Process EPA data
epa2020$Timestamp <- paste(epa2020$Date_GMT, epa2020$Hour_GMT)
epa2020$Timestamp <- paste(epa2020$Timestamp, ':00', sep = '')
col <- c('Lon', 'Lat', 'Timestamp', 'PM25')
epa <- epa2020 %>% 
  select(col)
epa <- epa[order(epa$Lon, epa$Lat, epa$Timestamp), ]

# Get 07/01/2020 - 07/07/2020 data
pa$Timestamp <- as.POSIXct(pa$Timestamp)
epa$Timestamp <- as.POSIXct(epa$Timestamp, format = '%Y-%m-%d %H:%M:%OS')
startTime <- as.POSIXct('2020-07-01 00:00:00')
endTime <- as.POSIXct('2020-07-01 23:59:59')
pa <- subset(pa, Timestamp <= endTime & Timestamp >= startTime)
epa <- subset(epa, Timestamp <= endTime & Timestamp >= startTime)

start <- as.numeric(startTime)
end <- as.numeric(as.POSIXct('2020-07-01 23:00:00'))
diff <- end - start

epa$indicator = 1
pa$indicator = 0

df_combine <- rbind(epa, pa)

# Construct NC prediction locations
s01   <- seq(-90,-65,0.1) # Lon
s02   <- seq(30,38,0.1) # Lat
s0    <- as.matrix(expand.grid(s01,s02))
innc <- map.where(map('state', 'north carolina', fill = TRUE), s0[,1], s0[,2])
s0    <- s0[!is.na(innc),]

# X0 is the prediction covariates
interim <- s0
interim <- dfColNorm(interim)
X0 <- locOpe(interim[, 1], interim[, 2])
X0 <- cbind(X0, 0)

n_timestamp <- diff / 3600 + 1
one_var <- vector(length = n_timestamp)

n.samples <- 25000
burn      <- 5000

# Model: Y = \beta_{0} + \beta_{1}*Lon + \beta_{2}*Lat + \beta_{3}*Lon^2 + 
# \beta_{4}*Lat^2 + \beta_{5}*LonLat + \beta_{6}I(epa)
for (i in 0:(n_timestamp - 1)) {
  
  current <- as.POSIXct(3600 * i, origin = startTime)
  
  
  subdf <- subset(df_combine, Timestamp == current)
  S <- data.matrix(subdf[, 1:2])
  maxd      <- max(dist(S))
  
  # Construct covarites
  interim <- S
  interim[, 1] <- (interim[, 1] - mean(interim[, 1])) / sd(interim[, 1])
  interim[, 2] <- (interim[, 2] - mean(interim[, 2])) / sd(interim[, 2])
  X <- locOpe(interim[, 1], interim[, 2])
  X <- cbind(X, subdf$indicator)
  
  
  
  # Kriging
  Y <- subdf$PM25
  
  starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(Y), "tau.sq"=0.5*var(Y))
  tuning    <- list("phi"=1, "sigma.sq"=0.05*var(Y), "tau.sq"=0.05*var(Y))
  amcmc     <- list("n.batch"=n.samples/100, "batch.length"=100, "accept.rate"=0.43)
  priors    <- list("beta.Norm"=list(rep(0,7), 100*diag(7)),
                    "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                    "sigma.sq.IG"=c(2, 1),
                    "tau.sq.IG"=c(2, 1))
  
  fit_spbayes  <- spLM(Y~X-1, coords=S,
                       starting=starting, tuning=tuning, 
                       priors=priors, cov.model="exponential",
                       amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
  
  pred <- spPredict(fit_spbayes, pred.coords=s0, pred.covars=X0, 
                    start=burn, thin=10, verbose=FALSE)
  
  pred <- pred$p.y.predictive.samples
  
  Yhat  <- apply(pred,1,mean)
  Ysd   <- apply(pred,1,sd)
  
  df <- data.frame(long=s0[,1],lat=s0[,2],Y=Ysd)
  # Samples near RTP area
  # Longitude range: -78.89, -78.52; Latitude range: 35.71, 35.92
  rtp_sample <- subset(df, df$long <=-78.52 & df$long >= -78.89 & df$lat >= 35.71
                       & df$lat <= 35.92)
  rtp_std <- mean(rtp_sample$Y)
  one_var[i + 1] = rtp_std
  # Plot standard deviation
  # Total number of pa monitor
  spa <- sum(subdf$indicator == 0)
  pred_var <- ggplot(df, aes(long, lat)) +
    geom_polygon(aes(x=long, y=lat)) +
    geom_raster(aes(fill = Y)) +
    scale_fill_gradientn(colours = viridis(10))+
    xlab(paste(current, 'pa:' ,spa))+ylab("")+labs(title="Posterior predictive standard deviation")+
    coord_fixed()
  plot(pred_var)
}

map <- map_data("county")
nc <- subset(map, region =="north carolina")
R <- ggplot(data=nc) + 
  geom_polygon(aes(x=long, y=lat, group = group)) +
  coord_fixed()
R

one_var

for (i in 1:7) {
  plot(one_var[(24*(i-1)+1):((24*i))], type = 'b')
}

plot(one_var, type = 'b')

