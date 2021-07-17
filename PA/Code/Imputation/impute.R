# Complete missing data in PA and EPA

library(tidyverse)
library(dplyr)
library(lutz)
library(lubridate)
rm(list=ls())

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Tools/myFunc.R')

pa2020 <- read.csv('Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv')
pa2020$Timestamp <- as.POSIXct(as.character(pa2020$Timestamp), format = "%Y-%m-%d %H:%M:%OS")
epa2020 <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv')
epa2020$Timestamp <- as.POSIXct(as.character(epa2020$Timestamp), format = "%Y-%m-%d %H:%M:%OS")

# Pivot PA data
pa <- pivot_wider(pa2020, names_from = c(Lon, Lat), values_from = PM25)
# Order by time
pa <- data.frame(pa[order(pa$Timestamp), ])

# Pivot EPA data
epa <- pivot_wider(epa2020, names_from = c(Lon, Lat), values_from = PM25)
# Order by time
epa <- data.frame(epa[order(epa$Timestamp), ])

## Complete data
# For the year of 2020, total number of timestamps shoud be 24 * 366 - 5= 8784 - 5 = 8779
# Remove monitors with less than 8000 readings leave us with 113 monitors
# As for the threshold (8000), this can be chosen using cross-validation later
dfna <- is.na(pa)
# A total of 142 monitors
comp <- apply(dfna, 2, sum)
newpa <- pa[, comp <= 8000]

# Complete date time
origin <- as.POSIXct('2020-01-01 00:00:00')
gap <- 3600 # 1 hr = 3600 sec
time2020 <- c()
for (i in 0:8778) {
  t <- as.POSIXct(i * gap, origin = origin)
  time2020 <- append(time2020, t)
}
timeframe <- data.frame(Timestamp = time2020)
newpa2 <- left_join(timeframe, newpa, by = 'Timestamp')
newpa <- data.frame(newpa2[order(newpa2$Timestamp), ])

# Now complete data on new pa
# Complete row 1 by averaging the row
colsum <- ncol(newpa)
row1 <- newpa[1, 2:colsum]
avg <- sum(row1, na.rm = TRUE) / sum(!is.na(row1))
row1[is.na(row1)] <- avg
newpa[1, 2:colsum] <- row1

# Use carry forward to complete the rest of the dataframe
imputeMissing <- function(a) {
  if (is.na(a[1]) == TRUE) {
    stop('Empty starting value')
  }
  for (i in 2:length(a)) {
    if (is.na(a[i]) == TRUE) {
      a[i] <- a[i - 1]
    }
  }
  return(a)
}

completePA <- apply(newpa, 2, imputeMissing)

# Pivot data frame back
completePA <- data.frame(completePA)
pa <- pivot_longer(completePA, cols = !Timestamp, names_to = 'coord', values_to = 'PM25')

temp <- pa %>% separate(coord, c('Lon', 'Lat'), sep = '_')
pa <- temp %>% separate(Lon, into = c(NA, 'Lon'), sep = 2)

pa$Timestamp <- as.POSIXct(pa$Timestamp)
pa$Lon <- -as.numeric(pa$Lon)
pa$Lat <- as.numeric(pa$Lat)
pa$PM25 <- as.numeric(pa$PM25)



## DO the same thing for epa data
# Complete date time
newepa2 <- left_join(timeframe, epa, by = 'Timestamp')
newepa <- data.frame(newpa2[order(newpa2$Timestamp), ])

colsum <- ncol(newepa)
row1 <- newepa[1, 2:colsum]
avg <- sum(row1, na.rm = TRUE) / sum(!is.na(row1))
row1[is.na(row1)] <- avg
newepa[1, 2:colsum] <- row1

completeEPA <- apply(newepa, 2, imputeMissing)

# Pivot data frame back
completeEPA <- data.frame(completeEPA)
epa <- pivot_longer(completeEPA, cols = !Timestamp, names_to = 'coord', values_to = 'PM25')

temp <- epa %>% separate(coord, c('Lon', 'Lat'), sep = '_')
epa <- temp %>% separate(Lon, into = c(NA, 'Lon'), sep = 2)

epa$Timestamp <- as.POSIXct(epa$Timestamp)
epa$Lon <- -as.numeric(epa$Lon)
epa$Lat <- as.numeric(epa$Lat)
epa$PM25 <- as.numeric(epa$PM25)

# Write
#write.csv(pa, 'PA_2020_Imputed.csv', row.names = FALSE)
#write.csv(epa, 'FRM_2020_Imputed.csv', row.names = FALSE)