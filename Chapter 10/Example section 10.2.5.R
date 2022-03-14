# This script shows the R code to create the contour
# shown in figure 10.4 part (a)
# for the example assuming GA distribution.

# Note: the llike function has 2 arguments.

library(gamlss)
data(aircond)

# log-likelihood function
llike <- function(mu, sigma) {
  sum(dGA(x=aircond, mu=mu, sigma=sigma, log=TRUE))
}

# likelihood function
like <- function(mu, sigma) {
  prod(dEXP(x=aircond, mu=mu, sigma=sigma))
}

# Vectorizing both functions
logL <- Vectorize(llike)
L    <- Vectorize(like)

# Creating one contour plot
mus    <- seq(from=40, to=120, by=0.1)
sigmas <- seq(from=0.6, to=1.6, by=0.05)
zz     <- outer(X=mus, Y=sigmas, FUN=logL)
contour(x=mus, y=sigmas, z=zz, nlevels=20,
        col=gray(0.3), lwd=2, lty='solid',
        xlab=expression(mu), ylab=expression(sigma))

