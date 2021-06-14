## Script used to format EPA and PA data
# Input: NC_PA_FRM_PM folder data
# Output: Formatted_PA_FRM data
rm(list=ls())
library(lubridate)
library(tidyverse)
library(dplyr)

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
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
epa2020$Timestamp <- paste(epa2020$Date, epa2020$Hour)
epa2020$Timestamp <- paste(epa2020$Timestamp, ':00', sep = '')
col <- c('Lon', 'Lat', 'Timestamp', 'PM25')
epa <- epa2020 %>% 
  select(col)
epa <- epa[order(epa$Lon, epa$Lat, epa$Timestamp), ]
pa$Timestamp <- as.POSIXct(pa$Timestamp)
epa$Timestamp <- as.POSIXct(epa$Timestamp, format = '%Y-%m-%d %H:%M:%OS')

# Convert Purple Air data from EST to GMT (5 hours difference)
pa$Timestamp <- with_tz(pa$Timestamp, 'UTC')
#write.csv(epa, 'FRM_2020_Hourly_Formatted.csv')
#write.csv(pa, 'PA_2020_Hourly_Formatted.csv')
