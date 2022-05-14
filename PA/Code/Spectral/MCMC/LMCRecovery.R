rm(list = ls())
library(fields) 
library(geoR)
library(truncnorm)
library(tidyverse)
library(mvtnorm)
library(gridExtra)
setwd("/Users/hongjianyang/Research/PAStudy/PA/Code/Spectral/MCMC/")
source('ExtraFunctions.R')
source('LMC_function.R')

#simulate data
source('simAllTS.R') # load your data here 

Y11 = Y1.real[, 1]
Y21 = Y2.real[, 1]
cor(Y11, Y21)
# names(out)=c('rangeU','rangeV','sigmaU','sigmaV','tau1','tau2','A','Y1.m','Y2.m','Y1p')

# Real data are Y1.real and Y2.real
# Coordinates are coords1 and coords2
iters = 3000
burn=1000
exit = LMC_fit(Y1=Y1.real,Y2=Y2.real, s1=coords1,s2=coords2,sp1=NULL,sp2=NULL,
               mean_range=0, sd_range=1, mean_var=0, sd_var=1, mean_rho=0,
               sd_rho=10, iters=iters, burn=burn, thin=1, update=10)

range1 = exit$rangeU
range2 = exit$rangeV
sigma1 = exit$sigmaU
sigma2 = exit$sigmaV
tau1 = exit$tau1
tau2 = exit$tau2
Al = exit$A

range11 = mean(range1)
range21 = range2[, 1]
sigma11 = sigma1[, , 1]
sigma21 = sigma2[, , 1]
sigma1mean = colMeans(t(sigma11))
sigma2mean = colMeans(t(sigma21))
tau11 = tau1[, 1]
tau21 = tau2[, 1]
A = Al[,, 1]

# Actual v.s. recovery
plot(1:3000, range1, 'l')
abline(rangeu,0,col=2,lwd=2)

plot(1:3000, range2, 'l')
abline(rangev,0,col=2,lwd=2)

plot(sigmau, sigma1mean)
abline(0,1,col=2,lwd=2)

plot(sigmav, sigma2mean)
abline(0,1,col=2,lwd=2)

plot(1:10, tau11)
abline(tau1,0,col=2,lwd=2)

plot(1:10, tau21)
abline(tau2,1,col=2,lwd=2)
