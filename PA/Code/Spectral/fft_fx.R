
fft_real <- function(dat,inverse=FALSE){
  if(!inverse){
    x  <- dat
    n  <- length(x)
    n2 <- floor(n/2)
    y  <- fft(x,inverse=FALSE)
    if(n%%2==0){
      X1     <- Re(y)[1:(n2+1)]
      X2     <- Im(y)[2:(n2)]
    }
    if(n%%2!=0){
      X1     <- Re(y)[1:(n2+1)]
      X2     <- Im(y)[2:(n2+1)]
    }
    out <- c(X1,X2)
  }
  if(inverse){
    X  <- dat
    n  <- length(X)
    n2 <- floor(n/2)
    if(n%%2==0){
      Y1    <- c(X[1:(n2+1)],X[n2:2])
      Y2    <- c(0,X[(n2+2):n],0,-X[n:(n2+2)])
    }
    if(n%%2!=0){
      Y1    <- c(X[1:(n2+1)],X[(n2+1):2])
      Y2    <- c(0,X[(n2+2):n],-X[n:(n2+2)])
    }
    y   <- complex(n, real = Y1, imaginary = Y2)
    out <- Re(fft(y/n,inverse=TRUE))
  }
return(out)}

# Test it

 n  <- 11
 Y  <- rnorm(n)
 Z  <- fft_real(Y,inverse=FALSE)
 Y2 <- fft_real(Z,inverse=TRUE)
 plot(Y2,Y);abline(0,1)

 n  <- 12
 Y  <- rnorm(n)
 Z  <- fft_real(Y,inverse=FALSE)
 Y2 <- fft_real(Z,inverse=TRUE)
 plot(Y2,Y);abline(0,1)

# Time it

 x <- rnorm(365)
 system.time(for(i in 1:10000){z<-fft_real(x)})

# Variance on spectral domain

 n <- 100000
 r <- 0.999
 Z <- rnorm(n)
 for(t in 2:n){Z[t] <- rnorm(1,r*Z[t-1],sqrt(1-r^2))}
 Y <- rnorm(n,Z,0.1)
 plot(Y-Z,main=var(Y-Z))

 y <- fft_real(Y)
 z <- fft_real(Z)
 plot((y-z)/sqrt(n/2),main=var((y-z)/sqrt(n/2)))
 plot(z)












