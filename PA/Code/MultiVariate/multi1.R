rm(list=ls())

library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis) 
library(lubridate)

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Tools/myFunc.R')




# type       := n-vector with type[i]=1 if Y[i] is an EPA measurement
#               and type[i]=2 if Y[i] is a PA measurement

# Read both EPA and PA data from Formatted_PA_FRM folder
# Output: PA, EPA
pa <- read.csv('Data/Formatted_PA_FRM/PA_2020_hourly_Formatted.csv')
epa <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv')
pa$Timestamp <- as.POSIXct(pa$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
epa$Timestamp <- as.POSIXct(epa$Timestamp, format = "%Y-%m-%d %H:%M:%OS")

# Train-test split for EPA stations
rownames(epa) <- NULL
set.seed(123)
epa_unique <- epa %>% count(Lon, Lat)
epa_unique$n <- 0

# Select four testing locations
epa_unique[2, 'n'] <- 1 # Near Asheville
epa_unique[6, 'n'] <- 1 # Near Charlotte
epa_unique[9, 'n'] <- 1 # Near Winston-Salem
epa_unique[16, 'n'] <- 1 # Near RDU and highway
epa_train_test <- left_join(epa, epa_unique, by = c('Lon', 'Lat'))
epa_train <- subset(epa_train_test, n == 0)
epa_test <- subset(epa_train_test, n == 1)


# Get a timestamp
start <- as.POSIXct('2020-05-07 22:00:00')
end <- as.POSIXct('2020-05-07 23:00:00')

startTime <- as.numeric(start)
endTime <- as.numeric(end)
diff <- endTime - startTime

n_timestamp <- diff / 3600 + 1

record <- c()
for (i in 0:(n_timestamp - 1)) {
  current <- as.POSIXct(3600 * i, origin = start)
  subepa_train <- subset(epa_train, Timestamp == current)
  subepa_test <- subset(epa_test, Timestamp == current)
  if (dim(subepa_test)[1] == 0) next
  subpa <- subset(pa, Timestamp == current)
  
  # subpa and subepa_train as training data, subepa_test as testing data
  # subpa first ,then subepa_train
  y <- rbind(data.frame(y = subpa$PM25), data.frame(y = subepa_train$PM25))$y
  s <- rbind(subpa[, c('Lon', 'Lat')], subepa_train[, c('Lon', 'Lat')])
  x <- s
  x <- dfColNorm(x)
  x <- locOpe(x[, 1], x[, 2])
  type1 <- data.frame(type = rep(1, nrow(subepa_train)))
  type2 <- data.frame(type = rep(2, nrow(subpa)))
  type <- rbind(type2, type1)
  
  fit <- LMC(y,x,s,type)
  
  
  
  # Do prediction at test locations
  theta_hat <- colMeans(fit$theta[1000:5000, ])
  d <- as.matrix(dist(s))
  t <- type
  sigma1 <- cov.LMC(d,type,theta_hat)$S
  sigma1_inv <- solve(sigma1)
  
  # Construct sigma0
  rho    <- theta_hat[1]
  range1 <- exp(theta_hat[2])
  range2 <- exp(theta_hat[3])
  tau1   <- exp(theta_hat[4])
  tau2   <- exp(theta_hat[5])
  sig1   <- exp(theta_hat[6])
  sig2   <- exp(theta_hat[7])
  
  s0 <- subepa_test[1, c('Lon', 'Lat')]
  s1 <- rbind(s0, s)
  # 1 * 90 distance matrix between s0 and s
  dist <- as.matrix(dist(s1))
  rownames(dist) <- NULL
  d <- t(as.matrix(dist[1,1:length(y)+1]))
  
  S      <- matrix(0,nrow = 1, ncol = length(y))
  S[1,type==1] <- sig1*exp(-d[1,type==1]/range1) +   
    sig2*rho*rho*exp(-d[1,type==1]/range2)
  S[1,type==2] <- sig2*rho*exp(-d[1,type==2]/range2)
  sigma0 <- S
  
  
  beta_hat <- colMeans(fit$beta[1000:5000, ])
  x0 <- dfColNorm(subepa_test[, c('Lon', 'Lat')])[1,]
  x0 <- locOpe(x0[1], x0[2])
  
  y_hat <- x0%*%beta_hat + sigma0 %*% sigma1_inv %*% (y - x %*% beta_hat)
  record <- append(record, y_hat)
}









