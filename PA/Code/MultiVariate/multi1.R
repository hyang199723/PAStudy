library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis) 
rm(list=ls())
setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Tools/myFunc.R')


## This is also a two-stage model
# In stage I, kriging on Purple Air only
# Get prediction readings in both epa_train
# In stage II, do a multivariate analysis and predict at epa_test

pa <- read.csv('Data/Formatted_PA_FRM/PA_2020_hourly_Formatted.csv')
epa <- read.csv('Data/Formatted_PA_FRM/FRM_2020_hourly_Formatted.csv')
pa$Timestamp <- as.POSIXct(pa$Timestamp)
epa$Timestamp <- as.POSIXct(epa$Timestamp)
## Global variables
# 80% of EPA stations were used during training
trainPercent = 0.8
# Number of covariates in prediction: Lon, Lat, lon^2, lat^2, lon*lat, intercept
p1 = 6

# Get 05/01/2020 - 05/01/2020 data
start <- as.POSIXct('2020-05-01 00:00:00')
end <- as.POSIXct('2020-05-01 00:00:00')
pa <- subset(pa, Timestamp >= start & Timestamp <= end)
epa <- subset(epa, Timestamp >= start & Timestamp <= end)

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

startTime <- as.numeric(start)
endTime <- as.numeric(end)
diff <- endTime - startTime

n_timestamp <- diff / 3600 + 1

n.samples <- 25000
burn      <- 5000

epa_train_sub <- epa_train
epa_test_sub <- epa_test
pa_sub <- pa


### Stage I: prediction with Purple Air data only
S <- data.matrix(pa_sub[, 1:2])
maxd      <- max(dist(S))
# Construct Purple Air covariates and observations
X <- S
X <- dfColNorm(X)
x.1 <- locOpe(X[, 1], X[, 2])
y.1 <- pa_sub$PM25
    
starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(Y), "tau.sq"=0.5*var(Y))
tuning    <- list("phi"=1, "sigma.sq"=0.05*var(Y), "tau.sq"=0.05*var(Y))
amcmc     <- list("n.batch"=n.samples/100, "batch.length"=100, "accept.rate"=0.43)
priors    <- list("beta.Norm"=list(rep(0,p1), 100*diag(p1)),
                      "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                      "sigma.sq.IG"=c(2, 1),
                      "tau.sq.IG"=c(2, 1))


stage1 <- spLM(y.1 ~ x.1 - 1, coords = S, starting = starting, tuning = tuning,
                   priors = priors, cov.model = "exponential", amcmc = amcmc, n.samples=n.samples,
                   verbose=FALSE)

# Predict at both epa_train locations
epa_train_s0 <- data.matrix(epa_train_sub[, 1:2])
epa_train_X0 <- epa_train_s0
epa_train_X0 <- dfColNorm(epa_train_X0)
epa_train_X0 <- locOpe(epa_train_X0[, 1], epa_train_X0[, 2])

pred_epa <- spPredict(stage1, pred.coords=epa_train_s0, pred.covars=epa_train_X0, 
                                       start=burn, thin=10, verbose=FALSE)

pred_epa_train <- pred_epa$p.y.predictive.samples
epa_train_hat <- apply(pred_epa_train,1,mean)


epa_test_s0 <- data.matrix(epa_test_sub[, 1:2])
epa_test_X0 <- epa_test_s0
epa_test_X0 <- dfColNorm(epa_test_X0)
epa_test_X0 <- locOpe(epa_test_X0[, 1], epa_test_X0[, 2])
    

#### Stage II: Multivariate analysis
# epa_train_hat: Krig

    # Construct EPA covariates and observations
    epa_s0 <- data.matrix(epa_train_sub[, 1:2])
    epa_X0 <- epa_s0
    x.2 <- dfColNorm(epa_X0)
    x.2 <- locOpe(x.2[, 1], x.2[, 2])
    y.2 <- epa_train_sub$PM25
    
    
    
    
    
    # Predict at testing locations
    lon <- epa_test_sub[, 1]
    lat <- epa_test_sub[, 2]
    test_s0 <- cbind(lon, lat)
    lon <- colNorm(lon)
    lat <- colNorm(lat)
    test_X0 <- locOpe(lon, lat)
    pred_test <- spPredict(fit_spbayes, pred.coords=test_s0, pred.covars=test_X0, 
                           start=burn, thin=10, verbose=FALSE)
    pred_test <- pred_test$p.y.predictive.samples
    epa_test_hat <- apply(pred_test,1,mean)


