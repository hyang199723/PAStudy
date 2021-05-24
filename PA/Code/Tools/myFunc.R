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