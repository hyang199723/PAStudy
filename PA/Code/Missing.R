library(tidyverse)
library(dplyr)
library(lutz)
library(lubridate)
rm(list=ls())

setwd("/Users/hongjianyang/Research/PAStudy/PA/")

pa2020 <- read.csv('Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv')
pa2020$Timestamp <- as.POSIXct(pa2020$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
epa2020 <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv')

# Pivot PA data
pa <- pivot_wider(pa2020, names_from = c(Lon, Lat), values_from = PM25)
pa <- 

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













