library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis)

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
source('Code/Tools/myFunc.R')
pa2020 <- read.csv('Data/NC_PA_FRM_PM/PA_2020_hourly_PM_NC.csv')
epa2020 <- read.csv('Data/NC_PA_FRM_PM/FRM_2020_hourly_PM_NC.csv')

## Global variables
# 80% of EPA stations were used during training
trainPercent = 0.8
# Number of covariates in stage 1: Lon, Lat, lon^2, lat^2, lon*lat, intercept
p1 = 6
# Number of covariates in stage 2: Lon, Lat, lon^2, lat^2, lon*lat, X(s), intercept
p2 = 7

col <- c('Lon', 'Lat', 'Timestamp', 'PM25_corrected')
pa <- pa2020 %>%
  select(col) %>%
  rename(PM25 = PM25_corrected) %>%
  group_by(Lon, Lat, Timestamp) %>%
  summarise(PM25 = mean(PM25))

pa <- pa[order(pa$Lon, pa$Lat, pa$Timestamp), ]
# Remove NAs
pa <- pa[complete.cases(pa), ]
pa$Timestamp <- as.POSIXct(pa$Timestamp)

# Process EPA data
epa2020$Timestamp <- paste(epa2020$Date_GMT, epa2020$Hour_GMT)
epa2020$Timestamp <- paste(epa2020$Timestamp, ':00', sep = '')
col <- c('Lon', 'Lat', 'Timestamp', 'PM25')
epa <- epa2020 %>% 
  select(col)
epa <- epa[order(epa$Lon, epa$Lat, epa$Timestamp), ]
epa$Timestamp <- as.POSIXct(epa$Timestamp, format = '%Y-%m-%d %H:%M:%OS')

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


    # Stage 1: Kriging on Purple Air data
    S <- data.matrix(pa_sub[, 1:2])
    maxd      <- max(dist(S))
    # Construct locations
    X <- S
    X <- dfColNorm(X)
    X <- locOpe(X[, 1], X[, 2])
    
    # Model: Y = \beta_{0} + \beta_{1}*Lon + \beta_{2}*Lat + \beta_{3}*Lon^2 + 
    # \beta_{4}*Lat^2 + \beta_{5}*LonLat
    Y <- pa_sub$PM25
    starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(Y), "tau.sq"=0.5*var(Y))
    tuning    <- list("phi"=1, "sigma.sq"=0.05*var(Y), "tau.sq"=0.05*var(Y))
    amcmc     <- list("n.batch"=n.samples/100, "batch.length"=100, "accept.rate"=0.43)
    priors    <- list("beta.Norm"=list(rep(0,p1), 100*diag(p1)),
                      "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                      "sigma.sq.IG"=c(2, 1),
                      "tau.sq.IG"=c(2, 1))
    
    fit_spbayes  <- spLM(Y~X-1, coords=S,
                         starting=starting, tuning=tuning, 
                         priors=priors, cov.model="exponential",
                         amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
    
    epa_s0 <- data.matrix(epa_train_sub[, 1:2])
    epa_X0 <- epa_s0
    epa_X0 <- dfColNorm(epa_X0)
    epa_X0 <- locOpe(epa_X0[, 1], epa_X0[, 2])
    
    # Predict at training locations
    pred_epa <- spPredict(fit_spbayes, pred.coords=epa_s0, pred.covars=epa_X0, 
                          start=burn, thin=10, verbose=FALSE)
    pred_epa <- pred_epa$p.y.predictive.samples
    epa_train_hat <- apply(pred_epa,1,mean)
    
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
    
    # Stage 2: Use Predicted values as covariates
    # Model: E(Y) = \beta_{0} + \beta_{1} * lon + \beta_{2} * lat + \beta_{3} * lon2
    #               + \beta_{4} * lat2 + \beta5 * lonlat + \beta6 * X(s)
    q = 2
    Y <- epa_train_sub$PM25
    maxd      <- max(dist(epa_s0))
    starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(Y), "tau.sq"=0.5*var(Y))
    tuning    <- list("phi"=1, "sigma.sq"=0.05*var(Y), "tau.sq"=0.05*var(Y))
    amcmc     <- list("n.batch"=n.samples/100, "batch.length"=100, "accept.rate"=0.43)
    priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
                   "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
    
    X.1 <- epa_X0
    X.2 <- epa_X0
    S <- epa_s0
    fit_spbayes_stage2  <- spMvLM(list(Y~X.1 - 1, epa_train_hat~X.2 - 1), coords=S,
                                starting=starting, tuning=tuning, 
                                priors=priors, cov.model="exponential",
                                amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
    


