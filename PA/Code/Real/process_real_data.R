                  #################################################
                  ################## Real data #################### 
                  #################################################

library(tidyverse)
library(dplyr)
library(maps)
library(spBayes)
library(viridis)
                
# Load data
pa2020 <- read.csv('./Real_data/Data/PA_2020_Hourly_Formatted.csv')
epa2020 <- read.csv('./Real_data/Data/FRM_2020_Hourly_Formatted.csv')

pa2020$Timestamp = as.POSIXct(pa2020$Timestamp)
epa2020$Timestamp = as.POSIXct(epa2020$Timestamp)

# set min and max date time
start.DT=max(min(pa2020$Timestamp),min(na.omit(epa2020$Timestamp)))
end.DT=as.POSIXct(c('2020-02-01 05:00:00 -03'))
#end.DT=min(max(pa2020$Timestamp),max(na.omit(epa2020$Timestamp)))

epa2020=epa2020 %>% filter(Timestamp>=start.DT&Timestamp<=end.DT)
pa2020=pa2020 %>% filter(Timestamp>=start.DT&Timestamp<=end.DT)

# Do these values make sense? 
# boxplot(epa2020$PM25)
# boxplot(pa2020$PM25)
# 
# range(epa2020$PM25,na.rm=TRUE)
# range(pa2020$PM25,na.rm=TRUE)

# create spatial id 

epa2020=epa2020%>% group_by(Lat,Lon) %>% mutate(id=cur_group_id())
n1=max(epa2020$id)
pa2020=pa2020%>% group_by(Lat,Lon) %>% mutate(id=cur_group_id()+n1)
n2=max(pa2020$id)-n1
ts=length(seq(start.DT, end.DT, by = "hour"))

pa2020 <- pa2020[order(pa2020$id,pa2020$Timestamp), ]
epa2020 <- epa2020[order(epa2020$id,epa2020$Timestamp), ]

## Complete columns with missing data data 
pa2020=pa2020 %>%
  complete(Timestamp = seq(start.DT, end.DT, by = "hour"))

epa2020=epa2020 %>%
  complete(Timestamp = seq(start.DT, end.DT, by = "hour"))

## refill Id sites and order
epa2020=epa2020%>% group_by(Lat,Lon) %>% mutate(id=cur_group_id())
pa2020=pa2020%>% group_by(Lat,Lon) %>% mutate(id=cur_group_id()+n1)
pa2020 <- pa2020[order(pa2020$id,pa2020$Timestamp), ]
epa2020 <- epa2020[order(epa2020$id,epa2020$Timestamp), ]

## Divide date and time
# pa2020$Date=as.Date(pa2020$Timestamp)
# pa2020$Time=format(pa2020$Timestamp,"%H:%M:%S")
# 
# epa2020$Date=as.Date(epa2020$Timestamp)
# epa2020$Time=format(epa2020$Timestamp,"%H:%M:%S")


## Join data
all.data=rbind(pa2020,epa2020)
all.data$Type=c(rep('PA',dim(pa2020)[1]),rep('EPA',dim(epa2020)[1]))


### Check NA proportions per site 
na.p=all.data %>% group_by(id,Type) %>% summarise(p.na=sum(is.na(PM25))/length(PM25))

ggplot(na.p)+
  geom_point(aes(x=id,y=p.na,col=Type))+theme_bw()+
  ggtitle('Proportions of NA values')

### Remove sites with more than 50% of NA values
bye.sites=which(na.p$p.na>=0.5)

all.dataF=all.data %>% filter(!id%in%bye.sites)

rm(all.data)

#check new proportions
na.pF=all.dataF %>% group_by(id,Type) %>% summarise(p.na=sum(is.na(PM25))/length(PM25))
ggplot(na.pF)+
  geom_point(aes(x=id,y=p.na,col=Type))+theme_bw()+
  ggtitle('Proportions of NA values')


sites.epa=all.dataF %>% filter(Type=='EPA') %>% summarise(ids=unique(id))
sites.epa=sites.epa$ids

sites.pa=all.dataF %>% filter(Type=='PA') %>% summarise(ids=unique(id))
sites.pa=sites.pa$ids

n1=length(sites.epa)
n2=length(sites.pa)

## Some plots 

# spatial correlation 
ggplot(all.dataF %>% filter(Timestamp==start.DT))+
  geom_point(aes(x=Lat,y=Lon,col=PM25))+facet_grid(~Type)+
  theme_bw()+coord_fixed(ratio = 1)+
  ggtitle('PM25 values')+scale_colour_gradient(low="#22FF00", high="#FF0000")

ggplot(all.dataF %>% filter(Timestamp==end.DT))+
  geom_point(aes(x=Lat,y=Lon,col=PM25))+facet_grid(~Type)+
  theme_bw()+coord_fixed(ratio = 1)+
  ggtitle('PM25 values')+scale_colour_gradient(low="#22FF00", high="#FF0000")

# Temporal correlation
### pa data
rs.pa=sample(sites.pa,size=10)
rs.epa=sample(sites.epa,size=10)

ggplot(all.dataF %>% filter(Type=='PA' & id%in%rs.pa))+geom_line(aes(x=Timestamp,y=PM25))+facet_wrap(~id)+theme_bw()
ggplot(all.dataF %>% filter(Type=='EPA' & id%in%rs.epa))+geom_line(aes(x=Timestamp,y=PM25))+facet_wrap(~id)+theme_bw()




