library(tidyverse)
library(dplyr)
library(lutz)
library(lubridate)
rm(list=ls())

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Tools/myFunc.R')
load('Data/Raw/PA2020_Raw1-50.Rdata')
a <- raw_data_full
load('Data/Raw/PA2020_Raw50-250.Rdata')
b <- raw_data_full
load('Data/Raw/PA2020_Raw250-400.Rdata')
c <- raw_data_full
load('Data/Raw/PA2020_Raw400-437.Rdata')
d <- raw_data_full

df <- rbind(a,b,c,d)
df1 <- df[order(df$CreatedAt),]
start <- as.POSIXct('2020-01-11 00:00:00')
end <- as.POSIXct('2020-01-11 23:00:00')

df2 <- subset(df1, (CreatedAt >= start & CreatedAt <= end))
df3 <- subset(pa2020, (Timestamp >= start & Timestamp <= end))

pa2020 <- read.csv('Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv')
pa2020$Timestamp <- as.POSIXct(as.character(pa2020$Timestamp))
epa2020 <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv')

# Pivot PA data
pa <- pivot_wider(pa2020, names_from = c(Lon, Lat), values_from = PM25)

# Get Missing Percentage
missing <- colSums(is.na(pa)) / nrow(pa)
hist(missing, main = '2020 Missing Hist')

cm.iterator <- getMonthTime()

for(i in 1:11) {
  sub_pa <- subset(pa2020, (Timestamp >= cm.iterator[i] & Timestamp < cm.iterator[i+1]))
  pa <- pivot_wider(sub_pa, names_from = c(Lon, Lat), values_from = PM25)
  sub_missing <- colSums(is.na(pa)) / nrow(pa)
  hist(sub_missing, main = paste('2020 Missing Hist', cm.iterator[i]))
}

# For December
sub_pa <- subset(pa2020, (Timestamp >= cm.iterator[12]))
pa <- pivot_wider(sub_pa, names_from = c(Lon, Lat), values_from = PM25)
sub_missing <- colSums(is.na(pa)) / nrow(pa)
hist(sub_missing, main = paste('2020 Missing Hist', cm.iterator[12]))


# Get number of monitors in each month
for(i in 1:11) {
  sub_pa <- subset(pa2020, (Timestamp >= cm.iterator[i] & Timestamp < cm.iterator[i+1]))
  pa <- pivot_wider(sub_pa, names_from = c(Lon, Lat), values_from = PM25)
  print(dim(pa))
}


##### Load Original PA data
old_pa <- read.csv('Data/Old/sensors.csv')












