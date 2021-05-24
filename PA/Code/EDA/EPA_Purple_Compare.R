library(mapproj)
library(usmap)
library(ggplot2)
library(viridis)
library(dplyr)
library(reshape2)
setwd("/Users/hongjianyang/Research/AirPollution/")
epa <- read.csv("Data/EPA.csv", header = TRUE)
p47479 <- read.csv("Data/Purple47479.csv", header = TRUE)
p83 <- read.csv("Data/Purple15583.csv", header = TRUE)
p87 <- read.csv("Data/Purple15587.csv", header = TRUE)
p91 <- read.csv("Data/Purple15591.csv", header = TRUE)

# Select columns
cols <- c('created_at', 'Temperature_F', 'Humidity_.', 'PM2.5_ATM_ug.m3')
p47479 <- subset(p47479, select = cols)
p83 <- subset(p83, select = cols)
p87 <- subset(p87, select = cols)
p91 <- subset(p91, select = cols)
# Truncate time
p47479$created_at <- trimws(vapply(strsplit(p47479$created_at," "), 
                                   `[`, 1, FUN.VALUE=character(1)))
p83$created_at <- trimws(vapply(strsplit(p83$created_at," "), 
                                   `[`, 1, FUN.VALUE=character(1)))
p87$created_at <- trimws(vapply(strsplit(p87$created_at," "), 
                                   `[`, 1, FUN.VALUE=character(1)))
p91$created_at <- trimws(vapply(strsplit(p91$created_at," "), 
                                   `[`, 1, FUN.VALUE=character(1)))


# Get daily average
mut_p47479 <- p47479 %>%
  rename(date = created_at, temperature = Temperature_F, humidity = Humidity_.,
         PM2.5 = PM2.5_ATM_ug.m3) %>%
  select(date, temperature, humidity, PM2.5) %>%
  group_by(date) %>%
  summarise(temperature = mean(temperature), humidity = mean(humidity), PM2.5 = mean(PM2.5))

mut_p83 <- p83 %>%
  rename(date = created_at, temperature = Temperature_F, humidity = Humidity_.,
         PM2.5 = PM2.5_ATM_ug.m3) %>%
  select(date, temperature, humidity, PM2.5) %>%
  group_by(date) %>%
  summarise(temperature = mean(temperature), humidity = mean(humidity), PM2.5 = mean(PM2.5))

mut_p87 <- p87 %>%
  rename(date = created_at, temperature = Temperature_F, humidity = Humidity_.,
         PM2.5 = PM2.5_ATM_ug.m3) %>%
  select(date, temperature, humidity, PM2.5) %>%
  group_by(date) %>%
  summarise(temperature = mean(temperature), humidity = mean(humidity), PM2.5 = mean(PM2.5))

# Pair EPA and Purple Air sensor
epa_pair <- epa[which(epa$Site.ID == 371830014), ]
data1 <- data.frame('date' = as.Date(mut_p47479$date), 'purple_2.5' = mut_p47479$PM2.5,
                   'epa_2.5' = epa_pair$Daily.Mean.PM2.5.Concentration)
data2 <- data.frame('date' = as.Date(mut_p83$date), 'purple_2.5' = mut_p83$PM2.5,
                    'epa_2.5' = epa_pair$Daily.Mean.PM2.5.Concentration)
data3 <- data.frame('date' = as.Date(mut_p87$date), 'purple_2.5' = mut_p87$PM2.5,
                    'epa_2.5' = epa_pair$Daily.Mean.PM2.5.Concentration)
df1 <- melt(data1, id.vars = 'date')
df2 <- melt(data2, id.vars = 'date')
df3 <- melt(data3, id.vars = 'date')
###### Generate Time Series Plot
# EPA v.s. Purple Air
ggplot(df1, aes(x=date, y=value)) +
  geom_line(aes(colour = variable)) +
  geom_point() +
  labs(title = 'PM2.5 Comparison near Raleigh (47479)', x = 'Date', 
       y = 'PM2.5 Concentration')
ggplot(df2, aes(x=date, y=value)) +
  geom_line(aes(colour = variable)) +
  geom_point() +
  labs(title = 'PM2.5 Comparison near Raleigh (83)', x = 'Date', 
       y = 'PM2.5 Concentration')
ggplot(df3, aes(x=date, y=value)) +
  geom_line(aes(colour = variable)) +
  geom_point() +
  labs(title = 'PM2.5 Comparison near Raleigh (87)', x = 'Date', 
       y = 'PM2.5 Concentration')
  
# Adjacent Purple Air sensors
data4 <- data.frame('date' = as.Date(mut_p87$date), '79_2.5' = mut_p47479$PM2.5,
                    '83_2.5' = mut_p83$PM2.5, '87_2.5' = mut_p87$PM2.5)
df4 <- melt(data4, id.vars = 'date')
ggplot(df4, aes(x=date, y=value)) +
  geom_line(aes(colour = variable)) +
  geom_point() +
  labs(title = 'PM2.5 Comparison between Purple Air sensors', x = 'Date', 
       y = 'PM2.5 Concentration')

###### Purple Air v.s. EPA 2020 monthly analysis
PA2020 <- read.csv("Data/PurpleAir_47479_2020.csv", header = TRUE)
EPA2020 <- read.csv("Data/EPA2020.csv", header = TRUE)

