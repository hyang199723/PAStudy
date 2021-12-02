# Compare speed of LMC functions
# Last update: 11/19/2021
rm(list=ls())
library(fields) 
library(geoR)
library(truncnorm)
library(tidyverse)
library(mvtnorm)
setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Spectral/ExtraFunctions.R')
source('Code/Spectral/LMC_function.R')

PA_data <- read.csv("Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv")
FRM_data <- read.csv("Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv")

# Convert timestamp
PA_data$Timestamp <- as.POSIXct(PA_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
FRM_data$Timestamp <- as.POSIXct(FRM_data$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
# No PA 
start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-01 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7

pa <- subset(PA_data, (Timestamp >= start) & (Timestamp <= end))
frm <- subset(FRM_data, (Timestamp >= start) & (Timestamp <= end))
# Get data to desired format
paTS <- pivot_wider(pa, names_from = Timestamp, values_from = PM25)
frmTS <- pivot_wider(frm, names_from = Timestamp, values_from = PM25)
# Record locations of PA and FRM stations
s1 <- as.matrix(frmTS[, 1:2])
s2 <- as.matrix(paTS[, 1:2])
# Get rid of the locations
paTS <- paTS[, -c(1:2)]
frmTS <- frmTS[, -c(1:2)]
Y1 = as.matrix(data.frame(frmTS))
colnames(Y1)=NULL
Y2 = as.matrix(data.frame(paTS))
colnames(Y2)=NULL

# Parameters
iters = 6000
thin = 1
# 1 Day speed
time1 <- proc.time()[3]
exit1 = LMC_fit(Y1,Y2, s1,s2, iters=iters, thin=thin)
time2 <- proc.time()[3]
elaps1 = time2 - time1
print(elaps1)




# 2 Day speed
start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-02 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7

pa <- subset(PA_data, (Timestamp >= start) & (Timestamp <= end))
frm <- subset(FRM_data, (Timestamp >= start) & (Timestamp <= end))
# Get data to desired format
paTS <- pivot_wider(pa, names_from = Timestamp, values_from = PM25)
frmTS <- pivot_wider(frm, names_from = Timestamp, values_from = PM25)
# Record locations of PA and FRM stations
s1 <- as.matrix(frmTS[, 1:2])
s2 <- as.matrix(paTS[, 1:2])
# Get rid of the locations
paTS <- paTS[, -c(1:2)]
frmTS <- frmTS[, -c(1:2)]
Y1 = as.matrix(data.frame(frmTS))
colnames(Y1)=NULL
Y2 = as.matrix(data.frame(paTS))
colnames(Y2)=NULL

time1 <- proc.time()[3]
exit1 = compact.LMC_fit(Y1,Y2, s1,s2, iters=iters, thin=thin)
time2 <- proc.time()[3]
elaps2 = time2 - time1
print(elaps2)



# 3 Day Speed
start = as.POSIXct('2020-03-01 05:00:00') 
end = as.POSIXct('2020-03-03 23:00:00') # 67 timstamps/spectrums Oct 2 FRM stations OCt 1 - 7

pa <- subset(PA_data, (Timestamp >= start) & (Timestamp <= end))
frm <- subset(FRM_data, (Timestamp >= start) & (Timestamp <= end))
# Get data to desired format
paTS <- pivot_wider(pa, names_from = Timestamp, values_from = PM25)
frmTS <- pivot_wider(frm, names_from = Timestamp, values_from = PM25)
# Record locations of PA and FRM stations
s1 <- as.matrix(frmTS[, 1:2])
s2 <- as.matrix(paTS[, 1:2])
# Get rid of the locations
paTS <- paTS[, -c(1:2)]
frmTS <- frmTS[, -c(1:2)]
Y1 = as.matrix(data.frame(frmTS))
colnames(Y1)=NULL
Y2 = as.matrix(data.frame(paTS))
colnames(Y2)=NULL

time1 <- proc.time()[3]
exit1 = compact.LMC_fit(Y1,Y2, s1,s2, iters=iters, thin=thin)
time2 <- proc.time()[3]
elaps3 = time2 - time1
print(elaps3)


