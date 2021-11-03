rm(list = ls())
# Plot EPA and PA sensors in North Carolina
# Last update: 11/01/2021
library(ggplot2)
library(maps)

setwd("/Users/hongjianyang/Research/PAStudy/PA/")

# Read in data
PA_data <- read.csv("Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv")
FRM_data <- read.csv("Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv")

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

  