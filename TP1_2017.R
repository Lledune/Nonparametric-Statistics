
## STAT 2150 ############# Class 1 (Oct. 26, 2017) #############

############################################################ clear all variables:
rm(list = ls())
library(stats)  #load packages (stats normally loaded by default)
library(datasets)

############################################################ Data generation

## Number of observations
n <- 100

## a) Simulation from a density 
# Z <- sort(runif(n, 0, 1)) #random number generator 
# X <- qnorm(Z) 
# plot(X,Z)

X <- rnorm(n, mean = 0, sd = 1)
plot(density(X))

## b) Random design m(x)
m <- function(x) (sin(2 * pi * x^3))^3
X <- sort(runif(n))
Y <- m(X) + 0.5 * rnorm(n, 0, 1)  #m(X) + iid-N(mean=0, sd=0.5) errors, added noise !
plot(X, Y, pch = "+")
rug(X) #points on the bottom
lines(X, m(X), col = 2)
# lines(X, Y)

## c) Fixed (equidistant) design : m(x)
X <- seq(from = 0, to = 1, length = n)
Y <- m(X) + 0.5 * rnorm(n, 0, 1)
plot(X, Y, pch = "+")
rug(X)
lines(X, m(X), col = 2)

## d) Data example: Cars (Random design)
X <- sort(cars$speed)
Y <- cars$dist[order(cars$speed)]
plot(X, Y, pch = "+")
rug(X)

############################################################ 
rm(list = ls())
## Working data
n <- 100
h <- 0.05  #normally need to determine h opt
X <- runif(n, 0, 1)  # observed data
x <- seq(from = min(X), to = max(X), length = 100)  # points on the domain 

############################################################ Different Kernels

Knorm <- function(u) dnorm(u)  #Gaussian kernel  
Kunif <- function(u) (abs(u) <= 1) * 0.5  #Uniform kernel
Kepan <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)  #Epanechnikov kernel
Ktria <- function(u) (1 - abs(u)) * (abs(u) <= 1)  #Triangle kernel

t <- seq(from = -2, to = 2, by = 0.005)
# Save the plot in the file pdf('D:/STAT2150/Plots/Plot1.pdf') # Pdf
# jpeg(file='D:/STAT2150/Plots/Plot1.jpeg') # JPEG
# png(file='D:/STAT2150/Plots/Plot1.png') # PNG

# pdf('D:/STAT2150/Plots/Plot1.pdf') # Pdf
## (pg. 21)
par(mfrow = c(2, 2)) 
plot(t, Knorm(t), type = "l", main = "Gaussian")
plot(t, Kunif(t), type = "l", main = "Uniform")
plot(t, Kepan(t), type = "l", main = "Epanechnikov")
plot(t, Ktria(t), type = "l", main = "Triangle")
dev.off() #shuts off graphs

##########################################################
## Probability density estimator
prDensEst <- function(x, X, h, K) mean(K((x - X)/h))/h

## Density kernel = 'gaussian', 'epanechnikov', 'rectangular' bw: the
## smoothing bandwidth The kernels are scaled such that bw is the
## standard deviation of the smoothing kernel
plot(density(X, from = min(X), to = max(X), bw = h, kernel = "gaussian"), 
    main = "Gaussian") #different bw
abline(h = 1, lty = 3)
rug(X)

# help(density) Density estimators
u2Knorm <- function(u) u^2 * dnorm(u)  # mu_2 Gaussian kernel  
u2Kunif <- function(u) u^2 * (abs(u) <= 1) * 0.5  # mu_2 Uniform kernel
u2Kepan <- function(u) u^2 * 0.75 * (1 - u^2) * (abs(u) <= 1)  # mu_2 Epanechnikov kernel

IntKnorm <- integrate(u2Knorm, -Inf, Inf)$value  # variance Gaussian kernel (= 1)
IntKunif <- integrate(u2Kunif, -Inf, Inf)$value  # variance Uniform kernel
IntKepan <- integrate(u2Kepan, -Inf, Inf)$value  # variance Epan kernel

h_norm <- h/IntKnorm^0.5 ## standardize bandwidth h to sd of Gaussian kernel
h_unif <- h/IntKunif^0.5 ## standardize bandwidth h to sd of Uniform kernel
h_epan <- h/IntKepan^0.5 ## standardize bandwidth h to sd of Epan kernal

fnorm <- sapply(x, function(x) prDensEst(x, X, h_norm, Knorm))
funif <- sapply(x, function(x) prDensEst(x, X, h_unif, Kunif))
fepan <- sapply(x, function(x) prDensEst(x, X, h_epan, Kepan))

## Gaussian kernel (1)
plot(density(X, n = 100, from = min(X), to = max(X), bw = h, kernel = "gaussian"), 
    main = "Kernel Density, Normal Kernel") ## built-in density estimator
lines(x, fnorm, col = 2) ## density estimator from scratch
abline(h = 1, lty = 3)
rug(X)  # the observations X

## Rectangular kernel (2)
plot(density(X, n = 100, from = min(X), to = max(X), bw = h, kernel = "rectangular"), 
    main = "Kernel Density, Uniform Kernel") ## built-in density estimator
lines(x, funif, col = 2) ## density estimator from scratch
abline(h = 1, lty = 3)
rug(X)  

## Epanechnikov kernel (3)
plot(density(X, n = 100, from = min(X), to = max(X), bw = h, kernel = "epanechnikov"), 
    main = "Kernel Density, Epanechnikov Kernel") ## built-in density estimator
lines(x, fepan, col = 2) ## density estimator from scratch
abline(h = 1, lty = 3)
rug(X)  

##########################################################
## Kernels and density at each point (pg. 18)
eachPoint <- function(X, h, K) 
{
    n <- length(X)
    f <- function(x) sapply(x, function(x) prDensEst(x, X, h, K))  #PR dens estim
    x <- seq(from = min(X), to = max(X), len = 100)
    ## plot points on x axis
    plot(X, rep(0, n), pch = "+", ylim = c(0, 1.25 * max(f(x))), main = "Kernel Density")
    ## draw a curve at each point
    for (Xi in X) {
        Ki <- K((x - Xi)/h)/(n * h)
        lines(x, Ki, type = "l", col = "blue")
    }
    # curve(f,add=T)
    lines(x, f(x))
}

# Figure of density estimators and kernels
# par(mfrow=c(2,2))
KnormAll <- eachPoint(X, h_norm, Knorm)  # Gaussian kernel
KepanAll <- eachPoint(X, h_epan, Kepan)  # Epanechnikov kernel
KunifAll <- eachPoint(X, h_unif, Kunif)  # Uniform kernel

############################################################
## Choice of bandwidth in practice
## Rule-of-thumb 'normal reference' bandwidth selector (pg. 38)

h_NR <- function(X, K) {
    RK <- integrate(function(u) K(u)^2, -Inf, Inf)$value ## integrated square kernel
    mu2 <- integrate(function(u) u^2 * K(u), -Inf, Inf)$value ## variance kernel (=1 for Gaussian)
    
    R <- quantile(X, 0.75) - quantile(X, 0.25) ## interquartile range
    sig <- min(sd(X), R/1.349)
    return(((8 * sqrt(pi) * RK)/(3 * mu2^2))^0.2 * sig * length(X)^(-0.2))
}

h1 <- h_NR(X, Knorm)      # rule-of-thumb 'normal reference' according to pg. 38
h3 <- bw.nrd(X)           # rule-of-thumb as implemented for "density" (bw.nrd = 1.06*sig*length(X)^(-0.2))
h4 <- bw.ucv(X, upper = 2*h1)     # unbiased cross-validation. implemented for "density" (pg. 40)
h6 <- bw.SJ(X)            # plug-in method of Sheather & Jones (1991) as implemented for "density" to select the
                          # bandwidth using pilot estimation of derivatives (pg. 39)

## Density with different bandwidths: 'nrd', 'ucv', 'SJ'
plot(range(X), c(0,1.5), type="n")
abline(h = 1, lty = 3)
lines(density(X, n=100, from=min(X), to=max(X), bw="nrd", kernel = "gaussian"), col=2)
lines(density(X, n=100, from=min(X), to=max(X), bw="ucv", kernel = "gaussian"), col=4)
lines(density(X, n=100, from=min(X), to=max(X), bw="SJ", kernel = "gaussian"), col=5)
legend("topleft", legend=c("nrd", "ucv", "SJ"),col = c(2, 4, 5), lty=c(1,1,1), cex=0.5)
rug(X)

##########################################################
## clear all variables:
rm(list = ls())

n <- 100
m <- function(x) sin((2 * pi * x^3))^3

X <- sort(runif(n, 0, 1)) # random design
Y <- m(X) + 0.5 * rnorm(n, mean = 0, sd = 1) # response variables
x <- seq(min(X), max(X), len = 100)  # points on the domain 
h <- 0.05

plot(X, Y, pch = "+")
lines(x, m(x), col = 2)
rug(X)

##########################################################
## NW regression estimator (pg. 48)
nwRegEst <- function(x, X, Y, h, K) sum(Y * K((x - X)/h))/sum(K((x - X)/h))

##########################################################
## Kernels cdf
cdfKnorm <- function(u) pnorm(u)  #Gaussian kernel cdf
cdfKunif <- function(u) (0.5 * (u + 1)) * (abs(u) <= 1) + (u > 1)  #Uniform kernel cdf
cdfKepan <- function(u) 0.25 * (2 + 3 * u - u^3) * (abs(u) <= 1) + (u > 1)  #Epanechnikov kernel cdf

##########################################################
## GM regression estimator (pg. 53)
gmRegEst <- function(x, X, Y, h, cdfK) {
    
    IntK <- NULL
    n <- length(X)
    S <- -Inf
    for (i in 1:(n - 1)) {
        S <- rbind(S, 0.5 * (X[i] + X[i + 1])) ## Compute s_i's (X already ordered)
    }
    S <- rbind(S, Inf)
    l <- (S - x)/h ## endpoints integrals
    
    for (i in 1:n) {
        IntK[i] <- cdfK(l[i + 1]) - cdfK(l[i]) ## Compute integrals
    }
    sum(IntK * Y)  # Y already ordered (because X is ordered)
}

#########################################################
h <- 0.02
Knorm <- function(u) dnorm(u)
NWnormEst <- sapply(x, function(x) nwRegEst(x, X, Y, h, Knorm))
GMnormEst <- sapply(x, function(x) gmRegEst(x, X, Y, h, cdfKnorm))

# par(mfrow=c(1,1))
# plot(X, Y, pch = "+")
plot(X, Y, pch=".")
lines(X, lm(Y ~ X)$fitted, col = 1)
lines(x, NWnormEst, type = "l", col = 4)
lines(x, GMnormEst, type = "l", col = 2)
lines(x, m(x), type = "l", col = 3)
legend("topleft", legend = c("LR", "NW", "GM", "True"), col = c(1, 4, 2, 3), 
       lty = c(1, 1, 1, 1), cex = 0.5)

#######################################
######## MONTE CARLO SIMULATION
#######################################

rm(list = ls())
K <- 1E4  # number of replicates
n <- 200  # number of observations

A <- matrix(nrow=K, ncol=2)

for (k in 1:K) {
    # set.seed(1)
    X <- cbind(rep(1, n), runif(n, 0, 1))
    Y <- X %*% c(2, 1.5) + 0.25 * rnorm(n, 0, 1)
    beta <- solve(t(X) %*% X) %*% t(X) %*% Y
    A[k, ] <- beta
}

plot(X[, 2], Y)  # last replicate
lines(X[, 2], mean(A[, 1]) + mean(A[, 2]) * X[, 2])  # overall fit
c(mean(A[, 1]), mean(A[, 2]))

(bias <- c(mean(A[, 1]) - 2, mean(A[, 2]) - 1.5))  # bias
(variance <- c(var(A[, 1]), var(A[, 2])))  # variance 
(ASE <- c(mean((A[, 1] - 2)^2), mean((A[, 2] - 1.5)^2)))  #MSE

hist(A[, 1], prob = T, col = "grey", breaks = 20)
abline(v = 2, lty = 2, lwd = 2, col = 2)
hist(A[, 2], prob = T, col = "grey", breaks = 20)
abline(v = 1.5, lty = 2, lwd = 2, col = 2)
