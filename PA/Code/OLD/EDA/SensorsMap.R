rm(list = ls())
# Plot EPA and PA sensors in North Carolina
# Last update: 11/05/2021
library(ggplot2)
library(maps)

setwd("/Users/hongjianyang/Research/PAStudy/PA/")

# Read in data
PA_data <- read.csv("Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv")
FRM_data <- read.csv("Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv")

PA_data$Timestamp <- as.POSIXct(PA_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
FRM_data$Timestamp <- as.POSIXct(FRM_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")

col_keep = c('Lon', 'Lat')
pa <- PA_data[, col_keep]
frm <- FRM_data[, col_keep]

# Get unique coordinates
pa <- unique(pa)
frm <- unique(frm)

# Plot
map <- map_data("county")
nc <- subset(map, region =="north carolina")
ggplot(data=nc) + 
  geom_polygon(aes(x=long, y=lat, group = group)) +
  geom_point(data = pa, aes(Lon, Lat), shape = 4, color = 'purple') +
  geom_point(data = frm, aes(Lon, Lat), shape = 16, color = 'red') + 
  labs(x = 'Lon', y = 'Lat', title = 'Distribution of Purple Air and FRM sensors in NC') +
  coord_fixed(ratio = 1)

# Time series plot of two stations
# PA station: -84.03269 35.06515 1/13
# FRM station: -83.44213 35.43477
PA_data$Lon = round(PA_data$Lon, 3)
PA_data$Lat = round(PA_data$Lat, 3)
FRM_data$Lon = round(FRM_data$Lon, 3)
FRM_data$Lat = round(FRM_data$Lat, 3)
station1 = subset(PA_data, (Lon == -84.033) & (Lat == 35.065))
station2 = subset(FRM_data, (Lon == -83.442) & (Lat == 35.435))

start = as.POSIXct('2020-01-23 00:00:00')
end = as.POSIXct('2020-01-23 23:00:00')

station1.sub = subset(station1, (Timestamp >= start) & (Timestamp <= end))
station2.sub = subset(station2, (Timestamp >= start) & (Timestamp <= end))
ggplot() +
  geom_line(data = station1.sub, aes(Timestamp, PM25), color = 'purple') +
  geom_line(data = station2.sub, aes(Timestamp, PM25), color = 'red') +
  labs(x = 'Time', y = 'PM2.5', title = 'Purple Air and EPA station reading comparison')

ggplot() +
  geom_line(data = station1, aes(Timestamp, PM25), color = 'purple') +
  labs(x = 'Time', y = 'PM2.5', title = 'One-year history readings of a Purple Air station')
  