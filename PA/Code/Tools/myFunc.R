# Given longitude and latitude, output second-order location matrix 
# including intercept
locOpe <- function(lon, lat) {
  lon2 <- lon * lon
  lat2 <- lat * lat
  lonlat <- lon * lat
  o <- cbind(1, lon, lat, lon2, lat2, lonlat)
  return(o)
}

# Normalize one column
colNorm <- function(col) {
  return ((col - mean(col, na.rm = TRUE)) / sd(col, na.rm = TRUE))
}
# Normalize a given data frame columnwise
dfColNorm <- function(df) {
  df <- apply(df, 2, colNorm)
  return(df)
}

# Get start and end time for each month in a calendar year
getMonthTime <- function() {
  jan.start <- as.POSIXct('2020-01-01 00:00:00')
  feb.start <- as.POSIXct('2020-02-01 00:00:00')
  mar.start <- as.POSIXct('2020-03-01 00:00:00')
  apr.start <- as.POSIXct('2020-04-01 00:00:00')
  may.start <- as.POSIXct('2020-05-01 00:00:00')
  jun.start <- as.POSIXct('2020-06-01 00:00:00')
  jul.start <- as.POSIXct('2020-07-01 00:00:00')
  aug.start <- as.POSIXct('2020-08-01 00:00:00')
  sep.start <- as.POSIXct('2020-09-01 00:00:00')
  oct.start <- as.POSIXct('2020-10-01 00:00:00')
  nov.start <- as.POSIXct('2020-11-01 00:00:00')
  dec.start <- as.POSIXct('2020-12-01 00:00:00')
  
  cm.iterator <- c(jan.start, feb.start, mar.start, apr.start, may.start, jun.start,
                   jul.start, aug.start, sep.start, oct.start, nov.start, dec.start)
  
  return(cm.iterator)
}
