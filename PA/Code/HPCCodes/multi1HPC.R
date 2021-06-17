rm(list=ls())

library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis) 
library(lubridate)

source('myFunc.R')




# type       := n-vector with type[i]=1 if Y[i] is an EPA measurement
#               and type[i]=2 if Y[i] is a PA measurement

# Read both EPA and PA data from Formatted_PA_FRM folder
# Output: PA, EPA
pa <- read.csv('PA_2020_Hourly_Formatted.csv')
epa <- read.csv('FRM_2020_Hourly_Formatted.csv')
pa$Timestamp <- as.POSIXct(pa$Timestamp, format = "%Y-%m-%d %H:%M:%OS")
epa$Timestamp <- as.POSIXct(epa$Timestamp, format = "%Y-%m-%d %H:%M:%OS")

# Normalize PA and EPA coordinates in one set, and split
full <- rbind(pa, epa)
full$LonNorm <- colNorm(full$Lon)
full$LatNorm <- colNorm(full$Lat)
split <- append(rep(1, nrow(pa)),rep(2, nrow(epa)))
pa <- subset(full, split == 1)
epa <- subset(full, split == 2)

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
start <- as.POSIXct('2020-05-01 00:00:00')
end <- as.POSIXct('2020-05-07 23:00:00')

startTime <- as.numeric(start)
endTime <- as.numeric(end)
diff <- endTime - startTime

n_timestamp <- diff / 3600 + 1

perform <- data.frame(Lon = epa_test$Lon, Lat = epa_test$Lat, 
                      Timestamp = epa_test$Timestamp, y = epa_test$PM25)
time <- c()
rho <- c()
logrange1 <- c()
logrange2 <- c()
logtau1 <- c()
logtau2 <- c()
logsig1 <- c()
logsig2 <- c()
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
  # Construct covariance matrix
  x <- rbind(subpa[, c('LonNorm', 'LatNorm')], subepa_train[, c('LonNorm', 'LatNorm')])
  x <- locOpe(x[, 1], x[, 2])
  type1 <- data.frame(type = rep(1, nrow(subepa_train)))
  type2 <- data.frame(type = rep(2, nrow(subpa)))
  type <- rbind(type2, type1)
  
  fit <- LMC(y,x,s,type)
  
  # Do prediction at test locations
  # Sample from posterior distribution in iteration 4k-5k
  record <- data.frame(ph = 1:nrow(subepa_test))
  for (j in 4001:5000) {
    theta_hat <- fit$theta[j, ]
    beta_hat <- fit$beta[j, ]
    s0 <- subepa_test[, c('Lon', 'Lat')]
    
    x0 <- subepa_test[, c('LonNorm', 'LatNorm')]
    x0 <- locOpe(x0[, 1], x0[, 2])
    y_hat <- LMCpredict(theta_hat, s, s0, x, x0, type, beta_hat, y)
    record <- cbind(record, y_hat)
  }
  record <- record[, 2:1001]
  
  # Construct time-theta matrix
  theta_hat <- apply(fit$theta[4001:5000, ], 2, mean)
  time <- append(time, current)
  rho <- append(rho, theta_hat[1])
  logrange1 <- append(logrange1, theta_hat[2])
  logrange2 <- append(logrange2, theta_hat[3])
  logtau1 <- append(logtau1, theta_hat[4])
  logtau2 <- append(logtau2, theta_hat[5])
  logsig1 <- append(logsig1, theta_hat[6])
  logsig2 <- append(logsig2, theta_hat[7])
  
  theta_temp <- apply(fit$theta[4001:5000, ], 2, mean)
  kriging_var <- LMCcovariance(theta_temp, s, s0, x, x0, type, beta_hat, y)
  
  # record is a n * 1000 matrix, where n is the number of testing locations
  temp <- data.frame(Lon = subepa_test$Lon, Lat = subepa_test$Lat, 
                     Timestamp = subepa_test$Timestamp, 
                     y_hat = apply(record, 1, mean), std = kriging_var^(1/2))
  
  perform <- left_join(perform, temp, by = c('Lon', 'Lat', 'Timestamp'))
  
  
}

theta_track <- data.frame(Timestamp = time, rho = rho, logrange1 = logrange1,
                          logrange2 = logrange2, logtau1 = logtau1, logtau2 = logtau2,
                          logsig1 = logsig1, logsig2 = logsig2)

write.csv(perform, 'perform.csv')
write.csv(theta_track, 'theta.csv')
