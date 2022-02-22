# Check imputation quality
rm(list = ls())
setwd('/Users/hongjianyang/Research/PAStudy/PA/')
library(raster)
library(rgdal)
roads.v <- readOGR("roads.shp")
plot(roads.v)
imputedata = read.csv("Data/Formatted_PA_FRM/EPA_Imputed_2020.csv")[,-1]
missingdata = read.csv("Data/Formatted_PA_FRM/missing_FRM.csv")[, -c(1,2,3)]

pa.missing <- read.csv("Data/Formatted_PA_FRM/missing_PA.csv")[, -c(1,2,3)]

stationone=3
stationtwo=4
missingd=t(missingdata[c(stationone,stationtwo),])
# extracting Row and Col positions of NA values
findNA=data.frame(which(is.na(missingd), arr.ind=TRUE))
findNA

X_stationone=findNA[which(findNA$col=='1'),]$row
X_stationone

X_stationtwo=findNA[which(findNA$col=='2'),]$row
X_stationtwo

imputed=t(imputedata[c(stationone,stationtwo),])
value1=imputed[X_stationone, 1]
value1
value2=imputed[X_stationtwo, 2]
value2

plot(missingd[, 1], type = 'l', main = 'Station 3',ylab = "PM2.5 readings")
points(x=X_stationone,y=value1,col="red")

plot(missingd[, 2], type = 'l', main = 'Station 4',ylab = "PM2.5 readings")
points(x=X_stationtwo,y=value2,col="red")


# Check stationary assumption of the original data
Y1 <- as.matrix(missingdata)
Y2 <- as.matrix(pa.missing)
Y1.missing <- is.na(Y1)
Y2.missing <- is.na(Y2)
# Mean imputation
beta1 <- mean(Y1, na.rm = T)
beta2 <- mean(Y2, na.rm = T)
Y1[Y1.missing] = beta1
Y2[Y2.missing] = beta2
# Look at residuals after de-trending
res1 = Y1 - beta1
res2 = Y2 - beta2

plot(res1[1,], type = 'l')
plot(res1[2,], type = 'l')
plot(res1[3,], type = 'l')

plot(res2[1,], type = 'l')
plot(res2[2,], type = 'l')
plot(res2[3,], type = 'l')

