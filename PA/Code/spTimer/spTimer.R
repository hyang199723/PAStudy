rm(list=ls())
setwd('/Users/hongjianyang/Research/PAStudy/PA/')

library(fields) 
library(glue)
library(viridis)
#library(geoR)
library(truncnorm)
library(tidyr)
library(mvtnorm)
library(ggplot2)
library(spTimer)
source('Code/Spectral/ExtraFunctions.R')
source('Code/Spectral/LMC_function.R')

#The data should be ordered first by the time and then by the sites specified by the coords below. 
# One can also supply coordi- nates through this argument, 
# where coordinate names should be "Latitude" and "Longitude".

# PM2.5 Data
pa <- read.csv('Data/Formatted_PA_FRM/PA_2020_Imputed.csv')
frm <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Imputed.csv')
pa$Timestamp <- as.POSIXct(levels(pa$Timestamp)[pa$Timestamp])
frm$Timestamp <- as.POSIXct(levels(frm$Timestamp)[frm$Timestamp])

# Start from a small subset: Use November data
startTime = as.POSIXct('2020-11-01')
endTime = as.POSIXct('2020-11-30 23:00:00')
pa.nov <- subset(pa, (Timestamp >= startTime) & (Timestamp <= endTime))
frm.nov <- subset(frm, (Timestamp >= startTime) & (Timestamp <= endTime))

# Mean-adjusting model: 
pa.nov$indicator <- 0
frm.nov$indicator <- 1

data <- rbind(pa.nov, frm.nov)
data <- data[with(data, order(Timestamp, Lon)), ]
data$Lon <- round(data$Lon, digits = 3)
data$Lat <- round(data$Lat, digits = 3)
rownames(data) <- NULL

data <- data[!duplicated(data), ]

testLon = -80.342 #sample(nrow(data), 1)
train <- subset(data, Lon != testLon)
test <- subset(data, Lon == testLon)

coords<-as.matrix(unique(cbind(train[,2:3])))
pred.coords<-as.matrix(unique(cbind(test[,2:3])))
# MCMC via Gibbs will provide output in *.txt format
# from C routine to avoide large data problem in R
set.seed(11)
post.gp.fitpred <- spT.Gibbs(formula = PM25 ~ Lon + Lat + indicator,
                             data=train, model="GP", coords=coords,
                             newcoords=pred.coords, newdata=test)
print(post.gp.fitpred)
summary(post.gp.fitpred)
coef(post.gp.fitpred)
plot(post.gp.fitpred)
names(post.gp.fitpred)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.gp.fitpred$prediction[,1]))
