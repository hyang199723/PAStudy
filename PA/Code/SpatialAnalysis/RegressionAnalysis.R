library(Stack) 
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
######## PM2.5 W.R.T temperature and humidity
df <- data.frame(pm2.5 = c(p47479$PM2.5_ATM_ug.m3, p83$PM2.5_ATM_ug.m3, p87$PM2.5_ATM_ug.m3,
                           p91$PM2.5_ATM_ug.m3),
                 tem = c(p47479$Temperature_F, p83$Temperature_F, p87$Temperature_F,
                         p91$Temperature_F),
                 humidity = c(p47479$Humidity_., p83$Humidity_., p87$Humidity_.,
                              p91$Humidity_.))

mls_both <- lm(df$pm2.5 ~ df$tem + df$humidity)
mls_tem <- lm(df$pm2.5 ~ df$tem)
mls_hum <- lm(df$pm2.5 ~ df$humidity)

mse1 <- sum((mls_both$fitted.values - df$pm2.5)^2) / (nrow(df) - 3)
mse2 <- sum((mls_tem$fitted.values - df$pm2.5)^2) / (nrow(df) - 2)
mse3 <- sum((mls_hum$fitted.values - df$pm2.5)^2) / (nrow(df) - 2)

########## EPA PM2.5 W.R.T Purple Air PM2.5
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

df1 <- data.frame(epa2.5 = epa_pair$Daily.Mean.PM2.5.Concentration,
                  pm79 = mut_p47479$PM2.5,
                  pm83 = mut_p83$PM2.5,
                  pm87 = mut_p87$PM2.5)

mls <- lm(df1$epa2.5 ~ df1$pm79 + df1$pm83 + df1$pm87)


