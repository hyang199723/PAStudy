library(tidyverse)
library(dplyr)

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Tools/myFunc.R')

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


