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
pa <- read.csv('Data/Formatted_PA_FRM/PA_2020_Hourly_Formatted.csv')
epa <- read.csv('Data/Formatted_PA_FRM/FRM_2020_Hourly_Formatted.csv')
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
start <- as.POSIXct('2020-05-01 05:00:00')
end <- as.POSIXct('2020-05-01 05:00:00')

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

#write.csv(perform, 'perform.csv')
#write.csv(theta_track, 'theta.csv')



#### Analyze perform.csv and theta.csv
perform <- read.csv('Outputs/MultiVariate/perform.csv')
theta <- read.csv('Outputs/MultiVariate/theta.csv')


# Clean raw perform.csv
if (F) {
perform[is.na(perform)] <- 0
perform$yhat <- 0
perform$std <- 0
for (i in 6:(ncol(perform) - 2)) {
  if (i %% 2 == 0) {
    perform$yhat <- perform$yhat + perform[, i]
  } else {
    perform$std <- perform$std + perform[, i]
  }
}
colKeep <- c('Lon', 'Lat', 'Timestamp', 'y', 'yhat', 'std')
perform <- select(perform, colKeep)
perform$Timestamp <- as.POSIXct(as.character(perform$Timestamp))
start <- as.POSIXct('2020-05-01 00:00:00')
end <- as.POSIXct('2020-05-07 23:00:00')

perform <- subset(perform, (Timestamp >= start) & (Timestamp <= end) )
}

# Total coverage probability
# 0.89
perform$lowerBound <- perform$yhat - qnorm(0.975) * perform$std
perform$upperBound <- perform$yhat + qnorm(0.975) * perform$std
perform$coverage <- as.numeric((perform$y >= perform$lowerBound) & (perform$y <= perform$upperBound))
sum(perform$coverage) / nrow(perform)
# Coverage probability by time
start <- as.POSIXct('2020-05-01 19:00:00')
end <- as.POSIXct('2020-05-01 20:00:00')

timeGroup <- perform %>% 
  group_by(Timestamp) %>%
  summarize(coverage = mean(coverage))

plot(x = as.POSIXct(timeGroup$Timestamp), y = timeGroup$coverage, xlab = 'Time', ylab = 'Coverage Probability')

perform$Lon = round(perform$Lon, 3)
# Coverage probability by location
ash <- subset(perform, Lon == -82.584)
sum(ash$coverage) / nrow(ash) #0.77

lotte <- subset(perform, Lon == -80.851)
sum(lotte$coverage) / nrow(lotte) #0.9

ws <- subset(perform, Lon == -80.342)
sum(ws$coverage) / nrow(ws) #0.98

rdu <- subset(perform, Lon == -78.820)
sum(rdu$coverage) / nrow(rdu) #0.89

# Calculate correlation
theta$var1 <- exp(theta$logsig1) + (theta$rho^2)*exp(theta$logsig2)+ exp(theta$logtau1)
theta$var2 <- exp(theta$logsig2)+ exp(theta$logtau2)
theta$corr <- theta$rho * exp(theta$logsig2) * ((sqrt(theta$var1) * sqrt(theta$var2))^(-1))

plot(theta$corr, type = 'b')
plot(theta$corr[1:24], type = 'b', xlab = 'day1')
plot(theta$corr[25:48], type = 'b', xlab = 'day2')
plot(theta$corr[49:72], type = 'b', xlab = 'day3')
plot(theta$corr[73:96], type = 'b', xlab = 'day4')
plot(theta$corr[97:120], type = 'b', xlab = 'day5')
plot(theta$corr[121:144], type = 'b', xlab = 'day6')
plot(theta$corr[145:168], type = 'b', xlab = 'day7')


plot(theta$corr[5:28], type = 'b', xlab = 'day1')
plot(theta$corr[29:52], type = 'b', xlab = 'day2')
plot(theta$corr[53:76], type = 'b', xlab = 'day3')
plot(theta$corr[77:100], type = 'b', xlab = 'day4')
plot(theta$corr[101:124], type = 'b', xlab = 'day5')
plot(theta$corr[121:144], type = 'b', xlab = 'day6')
plot(theta$corr[145:168], type = 'b', xlab = 'day7')



## Look at the posterior distribution of correlation
fitted_theta <- data.frame(fit$theta[4000:5000, ], check.names = T)
fitted_theta$var1 <- exp(fitted_theta$log.sig1) + (fitted_theta$rho^2)*exp(fitted_theta$log.sig2)+ exp(fitted_theta$log.tau1)
fitted_theta$var2 <- exp(fitted_theta$log.sig2)+ exp(fitted_theta$log.tau2)
fitted_theta$corr <- fitted_theta$rho * exp(fitted_theta$log.sig2) * ((sqrt(fitted_theta$var1) * sqrt(fitted_theta$var2))^(-1))
hist(fitted_theta$corr)



rm(list=ls())
n <- c(10, 100, 1000, 10000, 100000, 1e6)
variance <- c()
exp <- c()

fft_real <- function(dat,inverse=FALSE){
  if(!inverse){
    x  <- dat
    n  <- length(x)
    n2 <- floor(n/2)
    y  <- fft(x,inverse=FALSE)
    if(n%%2==0){
      X1     <- Re(y)[1:(n2+1)]
      X2     <- Im(y)[2:(n2)]
    }
    if(n%%2!=0){
      X1     <- Re(y)[1:(n2+1)]
      X2     <- Im(y)[2:(n2+1)]
    }
    out <- c(X1,X2)
  }
  if(inverse){
    X  <- dat
    n  <- length(X)
    n2 <- floor(n/2)
    if(n%%2==0){
      Y1    <- c(X[1:(n2+1)],X[n2:2])
      Y2    <- c(0,X[(n2+2):n],0,-X[n:(n2+2)])
    }
    if(n%%2!=0){
      Y1    <- c(X[1:(n2+1)],X[(n2+1):2])
      Y2    <- c(0,X[(n2+2):n],-X[n:(n2+2)])
    }
    y   <- complex(n, real = Y1, imaginary = Y2)
    out <- Re(fft(y/n,inverse=TRUE))
  }
  return(out)}


sim <- function(n) {
  delta <- rnorm(n, mean = 0, sd = 5)
  trans <- fft_real(delta)
  return(sd(trans)^2)
}

for (i in n) {
  v <- sim(i)
  variance <- append(variance, v)
  exp <- append(exp, i/2 * 25)
}

plot(x = n, y = log(variance), 'l')
lines(n, log(exp), col="green")


relD <- function(x,y) 2* abs(x - y) / abs(x + y)
















