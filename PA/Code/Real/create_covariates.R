##############################################################
#################### CREATE COVARIATES #######################
##############################################################

library(usmap)
library(ggplot2)
library(sp)
library(maps)
library(tidyverse)
library(rgeos)
library(sf)

# Load fire and smoke data 
load("./Real_data/HMS2021.RData")
names(fire_locs)=c('lon','lat','doy','hr')
fire_locs=data.frame(fire_locs)

# Load EPA and PA data
load("./Real_data/all_reail_init.RData") # update this 
all.dataF$Date=as.Date(all.dataF$Timestamp)
#all.dataF$Date=all.dataF$Date+366
all.dataF$Time=format(all.dataF$Timestamp,"%H")


ts=length(unique(all.dataF$Timestamp))
n1=length(sites.epa)  
n2=length(sites.pa)  

Ntot=length(all.dataF$Timestamp)
DTF=numeric(Ntot)
LS=numeric(Ntot)
MS=numeric(Ntot)
HS=numeric(Ntot)

# FIRE DATA
for (i in 1:Ntot)
{
  t=as.numeric(all.dataF[i,]$Time)
  d=all.dataF[i,]$Date
  d=as.numeric(d-as.Date("2021-01-01")+1)
  
  dataf=fire_locs %>% dplyr::filter(doy==d,hr==t)
  if (length(dataf$lon))
  {
    p <- tibble::tribble(~lon,~lat,all.dataF[i,]$Lon,all.dataF[i,]$Lat)
    coordinates(p)<-~lon+lat
    coordinates(dataf)<-~lon+lat
    DTF[i]=gDistance(p,dataf)  
  }
  else
  {
    DTF[i]=NA
  }

} 

# SMOKE DATA
sf::sf_use_s2(FALSE)

for (i in 1:Ntot)
{
  t=as.numeric(all.dataF[i,]$Time)
  d=all.dataF[i,]$Date
  d=as.numeric(d-as.Date("2021-01-01")+1)

  p <- as.data.frame(tibble::tribble(~lon,~lat,all.dataF[i,]$Lon,all.dataF[i,]$Lat))
  p_sf <- p %>% rowid_to_column() %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

  LS[i]=NA
  MS[i]=NA
  HS[i]=NA
  
  if(length(smoke_light[[d]]))
  {
    sl=st_sf(smoke_light[[d]])
    LS[i]=st_intersects(p_sf, sl,sparse = FALSE)
  }
  
  if(length((smoke_medium[[d]])))
  {
    sm=st_sf(smoke_medium[[d]])
    MS[i]=st_intersects(p_sf, sm,sparse = FALSE)  
  }
  
  if(length(smoke_heavy[[d]]))
  {
    sh=st_sf(smoke_heavy[[d]])
    HS[i]=st_intersects(p_sf, sh,sparse = FALSE)
  }
  
}


all.dataF$dtf=DTF
all.dataF$ls=LS
all.dataF$ms=MS
all.dataF$hs=HS




