# This code downloads data from this website
#      https://www.ospo.noaa.gov/Products/land/hms.html#data


rm(list=ls())
library(sf)

# Download HMS fire location data

filename<- function(m,d,y="2021"){
  dir <- "https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Fire_Points/Text/"
  f <- paste0(dir,y,"/",m,"/hms_fire",y,m,d,".txt")
return(f)}

mo  <- c("01","02","03","04","05","06","07","08","09","10","11","12")
da  <- c("01","02","03","04","05","06","07","08","09","10",
         "11","12","13","14","15","16","17","18","19","20",
         "21","22","23","24","25","26","27","28","29","30","31")

dat <- NULL
mos <- das <- NULL
doy <- 1
for(m in mo){for(d in da){
  fl   <- filename(m,d)
  temp <- try(as.matrix(read.csv(fl)), silent = TRUE,
              outFile = getOption("try.outFile", default = stderr()))
  if(is.matrix(temp)){
    print(c(m,d))
    mos[doy] <- m
    das[doy] <- d
    dat      <- rbind(dat,temp)
    doy      <- doy + 1
  }
}}

lon  <- as.numeric(dat[,1])
lat  <- as.numeric(dat[,2])
doy  <- as.numeric(dat[,3])-2021000
hr   <- floor(as.numeric(dat[,4])/100)
cali <- lon < -100 # Keep only these "near" California
lon  <- lon[cali]
lat  <- lat[cali]
doy  <- doy[cali]
hr   <-  hr[cali]

fire_locs <- cbind(lon,lat,doy,hr)

library(maps)
map("state")
points(lon,lat)



# Load smoke data
#   Plumes are classified as light, medium and heavy smoke
#   and these classes are downloaded separately


filename<- function(m,d,y="2021"){
  dir <- "https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/KML/"
  fl  <- paste0(dir,y,"/",m,"/hms_smoke",y,m,d,".kml")
return(fl)}

smoke_light <- list()
for(doy in 1:length(mos)){
  m    <- mos[doy]
  d    <- das[doy]
  fl   <- filename(m,d)
  temp <- try(st_read(fl,layer="Smoke (Light)"), silent = TRUE,
              outFile = getOption("try.outFile", default = stderr()))
  if(class(temp)[1]=="sf"){
    smoke_light[[doy]]  <- st_combine(temp)
    map("usa")
    plot(smoke_light[[doy]],add=TRUE,col=2)
    title(paste("month",m,"day",d))
  }
}


smoke_medium <- list()
for(doy in 1:length(mos)){
  m    <- mos[doy]
  d    <- das[doy]
  fl   <- filename(m,d)
  temp <- try(st_read(fl,layer="Smoke (Medium)"), silent = TRUE,
              outFile = getOption("try.outFile", default = stderr()))
  if(class(temp)[1]=="sf"){
    smoke_medium[[doy]]     <- st_combine(temp)
    map("usa")
    plot(smoke_medium[[doy]],add=TRUE,col=2)
    title(paste("month",m,"day",d))
  }
}


smoke_heavy <- list()
for(doy in 1:length(mos)){
  m    <- mos[doy]
  d    <- das[doy]
  fl   <- filename(m,d)
  temp <- try(st_read(fl,layer="Smoke (Heavy)"), silent = TRUE,
              outFile = getOption("try.outFile", default = stderr()))
  if(class(temp)[1]=="sf"){
    smoke_heavy[[doy]]     <- st_combine(temp)
    map("usa")
    plot(smoke_heavy[[doy]],add=TRUE,col=2)
    title(paste("month",m,"day",d))
  }
}

smoke_day   <- das
smoke_month <- mos
outfile     <- "/home/sofi/proyecto_doctoral/spatial_stat/PA_Data/Real_data/HMS2021.RData"
save(smoke_day,smoke_month,smoke_light,smoke_medium,smoke_heavy,fire_locs,file=outfile)





