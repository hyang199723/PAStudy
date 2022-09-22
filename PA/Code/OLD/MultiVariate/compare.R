rm(list=ls())

setwd("/Users/hongjianyang/Research/PAStudy/PA/")
stage <- read.csv('Outputs/TwoStage/perform.csv')
multi <- read.csv('Outputs/MultiVariate/perform.csv')
multi$bias <- multi$y - multi$yhat
multi$MSE <- multi$bias^2 + multi$std^2

# 95% convidence interval and converage
multi$lower <- multi$yhat - multi$std * 1.96 
multi$upper <- multi$yhat + multi$std * 1.96
multi$converage <- as.numeric((multi$lower <= multi$y) & (multi$upper >= multi$y))

# all stations from two-stage performance
stage1 <- subset(stage, Z == 100000)

# Overall coverage probability
sum(stage1$Cov) / nrow(stage1)   # 0.9412 for two-stage model
sum(multi$converage) / nrow(multi) # 0.892 for bivariate model

# coverage probability by station
stage1$Lon <- round(stage1$Lon, digits = 3)
multi$Lon <- round(multi$Lon, digits = 3)

ashe_stage <- subset(stage1, Lon == -82.584)
char_stage <- subset(stage1, Lon == -80.851)
ws_stage <- subset(stage1, Lon == -80.342)
rdu_stage <- subset(stage1, Lon == -78.820)

ashe_multi <- subset(multi, Lon == -82.584)
char_multi <- subset(multi, Lon == -80.851)
ws_multi <- subset(multi, Lon == -80.342)
rdu_multi <- subset(multi, Lon == -78.820)


sum(ashe_stage$Cov) / nrow(ashe_stage) # 0.945
sum(ashe_multi$converage) / nrow(ashe_multi)  #0.772

sum(char_stage$Cov) / nrow(char_stage) # 0.9125
sum(char_multi$converage) / nrow(char_multi)  #0.9

sum(ws_stage$Cov) / nrow(ws_stage) # 0.9994
sum(ws_multi$converage) / nrow(ws_multi)  #0.982

sum(rdu_stage$Cov) / nrow(rdu_stage) # 0.91
sum(rdu_multi$converage) / nrow(rdu_multi)  #0.89


# Compare overall MSE
mean(stage1$MSE) # 18.55
mean(multi$MSE) # 18.62

# MSE by location
mean(ashe_stage$MSE) # 28.337
mean(ashe_multi$MSE) # 17.17

mean(char_stage$MSE) # 18.76
mean(char_multi$MSE) # 21.23

mean(ws_stage$MSE) # 11.58
mean(ws_multi$MSE) # 14.47

mean(rdu_stage$MSE) # 16.74
mean(rdu_multi$MSE) # 21.74









 