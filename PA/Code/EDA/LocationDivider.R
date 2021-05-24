library(mapproj)
library(usmap)
library(ggplot2)
library(viridis)
setwd("/Users/hongjianyang/Research/AirPollution/")
sensors <- read.csv("Data/sensors.csv", header = TRUE)

# Get state for sensors. 11,700 US sensors
for (i in 1:nrow(sensors)){
  sensors[i, 'location'] <- map.where("state",x=sensors[i, 5], y=sensors[i, 4])
}
sensors$location <- gsub(':main', '', sensors$location)

# Get sensors in a specific state
state_sensors <- sensors[which(sensors$location == 'north carolina'), ]
# Export Data
write.csv(state_sensors, file = 'nc_sensors.csv')
############################################################################################
############################################################################################
############################################################################################
# US sensors
us_sensors <- sensors[which(!(is.na(sensors$location)) & sensors$pm2.5 < 100), ]

#######################
# Plot
ggplot(data = us_sensors, aes(longitude, latitude)) +
  borders('state', col = 'blue', fill = 'grey') + 
  geom_point(aes(color = pm2.5), size = 0.1) +
  scale_colour_gradientn(colours = viridis(20)) +
  coord_map(projection = "albers",lat0 = 39, lat1 = 45)

ggplot(data = state_sensors, aes(longitude, latitude)) +
  borders('state', col = 'blue', fill = 'grey') + 
  geom_point(aes(color = pm2.5), size = 1) +
  scale_colour_gradientn(colours = viridis(10)) +
  coord_map(projection = "albers",lat0 = 33, lat1 = 37,xlim=c(-85,-75),ylim=c(33,37))