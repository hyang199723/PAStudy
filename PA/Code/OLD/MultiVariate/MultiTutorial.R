library(spBayes)
rm(list=ls())
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu) 
  if(any(is.na(match(dim(V),p))))
    {stop("Dimension problem!")} 
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p))) 
}

set.seed(1)

# Generate some data
n <- 25 # number of locations
q <- 2 # number of outcomes at each location
nltr <- q*(q+1)/2 # number of triangular elements in the cross-covariance matrix
# matrix of observation coordinates for spMvLM()
coords <- cbind(runif(n,0,1), runif(n,0,1))
# Parameters for the bivariate spatial random effects
theta <- rep(3/0.5,q)
# Lower-triangular matrix for cross-covariance of Gaussian process
A <- matrix(0,q,q) 
A[lower.tri(A,TRUE)] <- c(1,-1,0.25)
K <- A%*%t(A)


# dispersion matrix
Psi <- diag(0,q)
# calculate spatial covariance matrix
C <- mkSpCov(coords, K, Psi, theta, cov.model="exponential") # Gaussian spatial process
w <- rmvn(1, rep(0,nrow(C)), C)
# w.1 and w.2 are every other element in w
w.1 <- w[seq(1,length(w),q)] 
w.2 <- w[seq(2,length(w),q)]
# Covariate portion of the mean
x.1 <- cbind(1, rnorm(n))
x.2 <- cbind(1, rnorm(n))
# create a multivariate design matrix given q univariate design matrices 
x <- mkMvX(list(x.1, x.2))
# parameters
B.1 <- c(1,-1) 
B.2 <- c(-1,1)
B <- c(B.1, B.2)
# dispersion
Psi <- diag(c(0.1, 0.5)) # outcomes
y <- rnorm(n*q, x%*%B+w, diag(n)%x%Psi)
# outcomes based on every other element in y
y.1 <- y[seq(1,length(y),q)] 
y.2 <- y[seq(2,length(y),q)]



# Call spMvLM
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)] # number of MCMC iterations for spMvLM()
n.samples <- 1000
# tags with starting values for parameters for spMvLM()
starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
# tags with Metropolis sampler variances for spMvLM()
tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q)) # tags with prior values
priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
               "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
# call spMvLM() function
m.1 <- spMvLM(list(y.1~x.1-1, y.2~x.2-1),
              coords=coords, starting=starting, tuning=tuning, priors=priors, n.samples=n.samples, cov.model="exponential", n.report=100)
burn.in <- 0.75*n.samples
m.1 <- spRecover(m.1, start=burn.in)
