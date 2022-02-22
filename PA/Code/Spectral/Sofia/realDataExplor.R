##
# This is Sofia's codes to look at real data
# It is mainly EDA: data cleaning and visualization


##################################################################################################################
###################################### Real Data exploration  ##################################################
##################################################################################################################

library(geoR)
library(ggplot2)
library(tidyverse)
library(sp)


###
## 1. Load data
###

# load purple air data
PA_data <- read.csv("./Real_data/PA_2020_Hourly_Formatted.csv")
# load Environmental Protection agency data
FRM_data <- read.csv("./Real_data/FRM_2020_Hourly_Formatted.csv")

###
## 2. Tidy data:
###

# 2.1 Converting to date time
PA_data$DateTime=as.POSIXct(PA_data$Timestamp,tz=Sys.timezone())
FRM_data$DateTime=as.POSIXct(FRM_data$Timestamp,tz=Sys.timezone())

# 2.2 extracting time
PA_data$Time <- format(as.POSIXct(PA_data$Timestamp),format = "%H:%M:%S")
PA_data$Time <- as.POSIXct(PA_data$Time,format="%H:%M:%S")

FRM_data$Time <- format(as.POSIXct(FRM_data$Timestamp),format = "%H:%M:%S")
FRM_data$Time <- as.POSIXct(FRM_data$Time,format="%H:%M:%S")

# 2.3 extracting Date
PA_data$Date <- as.Date(PA_data$Timestamp)
FRM_data$Date <- as.Date(FRM_data$Timestamp)

# 2.4  Group data by site 
nsitesPA=length(unique(PA_data$Lon+PA_data$Lat))
nsitesFRM=length(unique(FRM_data$Lon+FRM_data$Lat))
PA_data=PA_data %>% group_by(Lon+Lat) %>%  mutate("idSite" = as.factor(Lon+Lat))
FRM_data= FRM_data %>% group_by(Lon+Lat) %>%  mutate("idSite" = as.factor(Lon+Lat))
levels(PA_data$idSite)=seq(1:nsitesPA)
levels(FRM_data$idSite)=seq(1:nsitesFRM)

# number of observations per site
n.by.sitePA=PA_data %>% group_by(idSite) %>% summarise(n=n())
n.by.siteFRM=FRM_data %>% group_by(idSite) %>% summarise(n=n())

# It seems some locations have only one observation...
# remove them
whichId=PA_data %>% group_by(idSite) %>% summarise(n=n()) %>% subset(n!=1) %>% select(idSite)
PA_filter=PA_data %>% filter(idSite%in%whichId$idSite)

# Calculate mean point
mean_lonPA=(min(PA_filter$Lon)+max(PA_filter$Lon))/2
mean_latPA=(min(PA_filter$Lat)+max(PA_filter$Lat))/2
mean_lonFRM=(min(FRM_data$Lon)+max(FRM_data$Lon))/2
mean_latFRM=(min(FRM_data$Lat)+max(FRM_data$Lat))/2

Dist_to_center=spDistsN1(pts=cbind(PA_filter$Lon,PA_filter$Lat),pt=c(mean_lonPA,mean_latPA),longlat=TRUE)
PA_filter$dtc=Dist_to_center
Dist_to_center_FRM=spDistsN1(pts=cbind(FRM_data$Lon,FRM_data$Lat),pt=c(mean_lonFRM,mean_latFRM),longlat=TRUE)
FRM_data$dtc=Dist_to_center_FRM



###
## 3. Prelimirary plots
###

# 3.1 spatail distribution 
dtime="2020-01-12 18:00:00"
ggplot(PA_data %>% filter(Timestamp==dtime))+
  geom_point(aes(Lon, Lat, color = PM25), size = 2) +
  coord_fixed(ratio = 1) +
  scale_color_gradient(low = "blue", high = "orange") +
  borders("county", "north carolina")+
  theme_bw()

ggplot(FRM_data %>% filter(Timestamp==dtime))+
  geom_point(aes(Lon, Lat, color = PM25), size = 2) +
  coord_fixed(ratio = 1) +
  scale_color_gradient(low = "blue", high = "orange") +
  borders("county", "north carolina")+
  theme_bw()

# 3.2 number of observations per site 
plot(n.by.sitePA)
plot(n.by.siteFRM)
boxplot(n.by.sitePA$n)
boxplot(n.by.siteFRM$n)


# 3.3: Time series by site  
# PA data
ggplot(PA_filter %>% filter(idSite%in%seq(4,20)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)
ggplot(PA_filter %>% filter(idSite%in%seq(21,40)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)
ggplot(PA_filter %>% filter(idSite%in%seq(41,60)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)
ggplot(PA_filter %>% filter(idSite%in%seq(61,80)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)
ggplot(PA_filter %>% filter(idSite%in%seq(81,100)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)
ggplot(PA_filter %>% filter(idSite%in%seq(101,120)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)
ggplot(PA_filter %>% filter(idSite%in%seq(121,141)))+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)

#FRM data
ggplot(FRM_data)+geom_line(aes(x=DateTime,y=PM25,color=dtc))+theme_bw()+facet_wrap(~idSite)

