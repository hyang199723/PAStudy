library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis)


## Global variables
# 80% of EPA stations were used during training
trainPercent = 0.8
# Number of covariates in stage 1: Lon, Lat, lon^2, lat^2, lon*lat, intercept
p1 = 6
# Number of covariates in stage 2: Lon, Lat, lon^2, lat^2, lon*lat, X(s), intercept
p2 = 7

setwd("/Users/hongjianyang/Research/PA/")
source('Code/Tools/myFunc.R')
pa2020 <- read.csv('Data/NC_PA_FRM_PM/PA_2020_hourly_PM_NC.csv')
epa2020 <- read.csv('Data/NC_PA_FRM_PM/FRM_2020_hourly_PM_NC.csv')
# Subset data
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



# Get 05/05/2020 data

start <- as.POSIXct('2020-05-05 00:00:00')
end <- as.POSIXct('2020-05-05 23:00:00')
pa <- subset(pa, Timestamp >= start & Timestamp <= end)
epa <- subset(epa, Timestamp >= start & Timestamp <= end)

### Two Stage ALgorithm
## In stage 1, kriging on Purple Air only.
# X: covariate matrix, p = 6, purple air locations *
# S: Purple Air locations *
# epa_X0: covariates when predicting Pm25 at EPA stations (p=6) *
# epa_s0: locations of epa stations *
# EPA stations in stage1 would be 80% of all data

## In stage 2, kriging on 80% epa stations
# X: covariate matrix, p = 7, EPA locations and predicted values from stage 1 *
#     (a column bind of epa_X0 and predicted value from 1)
# S: EPA locations *
#     (the same as epa_s0 from 1)
# test_X0: covariate matrix for 20% EPA stations *
# test_s0: 20% EPA stations to be predicted on *

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
performTwo <- data.frame()
for (i in 0:(n_timestamp - 1)) {
  current <- as.POSIXct(3600 * i, origin = start)
  epa_train_sub <- subset(epa_train, Timestamp == current)
  epa_test_sub <- subset(epa_test, Timestamp == current)
  if (dim(epa_test_sub)[1] == 0) next
  pa_sub <- subset(pa, Timestamp == current)
  
  z <- c(1000, 5,4,3,2,1)
  avg = mean(pa_sub$PM25)
  deviation = sd(pa_sub$PM25)
  epa_avg <- mean(epa_train_sub$PM25)
  epa_std <- sd(epa_train_sub$PM25)
  print(i/n_timestamp)
  for (j in z) {
    # Exclude extreme values
    pa_sub_z <- subset(pa_sub, (PM25 <= avg + j*deviation) 
                       & (PM25 >= avg - j*deviation))
    epa_sub_z <- subset(epa_train_sub, (PM25 <= epa_avg + j*epa_std) 
                        & (PM25 >= epa_avg - j*epa_std))
    # Stage 1: Kriging on Purple Air data
    S <- data.matrix(pa_sub_z[, 1:2])
    maxd      <- max(dist(S))
    # Construct locations
    X <- S
    X <- dfColNorm(X)
    X <- locOpe(X[, 1], X[, 2])
    
    # Model: Y = \beta_{0} + \beta_{1}*Lon + \beta_{2}*Lat + \beta_{3}*Lon^2 + 
    # \beta_{4}*Lat^2 + \beta_{5}*LonLat
    Y <- pa_sub_z$PM25
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
    
    epa_s0 <- data.matrix(epa_sub_z[, 1:2])
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
    Y <- epa_sub_z$PM25
    maxd      <- max(dist(epa_s0))
    starting  <- list("phi"=1/(0.05*maxd), "sigma.sq"=0.5*var(Y), "tau.sq"=0.5*var(Y))
    tuning    <- list("phi"=1, "sigma.sq"=0.05*var(Y), "tau.sq"=0.05*var(Y))
    amcmc     <- list("n.batch"=n.samples/100, "batch.length"=100, "accept.rate"=0.43)
    priors    <- list("beta.Norm"=list(rep(0,p2), 100*diag(p2)),
                      "phi.Unif"=c(1/(2*maxd), 1/(0.01*maxd)),
                      "sigma.sq.IG"=c(2, 1),
                      "tau.sq.IG"=c(2, 1))
    
    X <- cbind(epa_X0, epa_train_hat)
    S <- epa_s0
    fit_spbayes_stage2  <- spLM(Y~X-1, coords=S,
                                starting=starting, tuning=tuning, 
                                priors=priors, cov.model="exponential",
                                amcmc=amcmc, n.samples=n.samples,verbose=FALSE)
    
    test_X0 <- cbind(test_X0, epa_test_hat)
    pred_stage2 <- spPredict(fit_spbayes_stage2, pred.coords=test_s0, pred.covars=test_X0, 
                             start=burn, thin=10, verbose=FALSE)
    
    pred_stage2 <- pred_stage2$p.y.predictive.samples
    Yhat <- apply(pred_stage2, 1, mean)
    std <- apply(pred_stage2, 1, sd)
    Y <- epa_test_sub$PM25
    
    MSE <- (Y - Yhat)^2 + std^2
    subRec <- select(epa_test_sub, c('Lon', 'Lat', 'Timestamp'))
    subRec$bias <- Y - Yhat
    subRec$std <- std
    subRec$MSE <- MSE
    subRec$Cov <- as.numeric((Y >= Yhat-1.96*std) & (Y <= Yhat+1.96*std))
    subRec$Z <- j
    performTwo <- rbind(performTwo, subRec)
  }
}

#beta <- spRecover(fit_spbayes_stage2, start=burn, verbose=FALSE)$p.beta.recover.samples
#apply(beta, 2, mean)
performTwo$Lon <- round(performTwo$Lon, 3)

# -82.58440, -80.85147, -80.34200, -78.81970
# Check performance near Asheville
ashe <- subset(performTwo, Lon == -82.584)
char <- subset(performTwo, Lon == -80.851)
wins <- subset(performTwo, Lon == -80.342)
rdu <- subset(performTwo, Lon == -78.820)

ggplot(data = ashe, aes(x=Timestamp, y=MSE, group=Z, color=Z)) +
  geom_line() +
  scale_color_viridis()

ggplot(data = ashe, aes(x=Timestamp, y=bias, group=Z, color=Z)) +
  geom_line() +
  scale_color_viridis()

ggplot(data = ashe, aes(x=Timestamp, y=std, group=Z, color=Z)) +
  geom_line() +
  scale_color_viridis()

ashe1000 <- subset(ashe, Z == 1000)
ashe1 <- subset(ashe, Z == 1)
ashe3 <- subset(ashe, Z == 3)
plot(x = ashe1000$Timestamp, y = ashe1000$MSE, 'b')
plot(x = ashe1$Timestamp, y = ashe1$MSE,  'b')
plot(x = ashe3$Timestamp, y = ashe3$MSE,  'b')
plot(x = ashe3$Timestamp, y = ashe3$bias,  'b')
plot(x = ashe3$Timestamp, y = ashe3$std,  'b')

# Coverage Probability
asheCov <- ashe[ashe$Z == 1000, ]
sum(asheCov$Cov) / nrow(asheCov) # 0.8695652, z = 1000
asheCov <- ashe[ashe$Z == 5, ]
sum(asheCov$Cov) / nrow(asheCov) # 0.8695652, z = 5
asheCov <- ashe[ashe$Z == 4, ]
sum(asheCov$Cov) / nrow(asheCov) # 0.8695652, z = 4
asheCov <- ashe[ashe$Z == 3, ]
sum(asheCov$Cov) / nrow(asheCov) # 0.9130435, z = 3
asheCov <- ashe[ashe$Z == 2, ]
sum(asheCov$Cov) / nrow(asheCov) # 0.8695652, z = 2
asheCov <- ashe[ashe$Z == 1, ]
sum(asheCov$Cov) / nrow(asheCov) # 0.6956522, z = 1


charCov <- char[char$Z == 3, ]
sum(charCov$Cov) / nrow(charCov) # 0.9130435, z = 3
