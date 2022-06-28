# This script shows the R code to use checkMomentsSK function
# with simulated data in which we know the distribution

# Example with NO ---------------------------------------------------------

x <- rNO(n=1000, mu=170, sigma=10)
momentSK(x)

# Using the random sample
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Using the best fit from fitDist
m1 <- fitDist(x, type="realline")
m1$fits
checkMomentSK(m1)

# Using mod from gamlss function
# assuming the correct distribution
mod <- gamlss(x ~ 1, family=NO)
checkMomentSK(mod)

# Using mod from gamlss function
# assuming a wrong distribution
mod <- gamlss(x ~ 1, family=WEI)
checkMomentSK(mod)

# Note: when assuming the correct model we can 
# obtain a plot in which "mod" is located
# inside ellipse


# Example with GA ---------------------------------------------------------

x <- rGA(n=1000, mu=17, sigma=0.5)
momentSK(x)

# Using the random sample
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Using the best fit from fitDist
m1 <- fitDist(x, type="realplus")
m1$fits
checkMomentSK(m1)

# Using mod from gamlss function
# assuming the correct distribution
mod <- gamlss(x ~ 1, family=GA)
checkMomentSK(mod)

# Using mod from gamlss function
# assuming a wrong distribution
mod <- gamlss(x ~ 1, family=WEI)
checkMomentSK(mod)

# Note: when assuming the correct model we can 
# obtain a plot in which "mod" is located
# inside ellipse


# Example with PO ---------------------------------------------------------

x <- rPO(n=1000, mu=170)
momentSK(x)

# Using the random sample
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Using the best fit from fitDist
m1 <- fitDist(x, type="counts")
m1$fits
checkMomentSK(m1)

# Using mod from gamlss function
# assuming the correct distribution
mod <- gamlss(x ~ 1, family=PO)
checkMomentSK(mod)

# Using mod from gamlss function
# assuming a wrong distribution
mod <- gamlss(x ~ 1, family=GEOM)
checkMomentSK(mod)

# Note: when assuming the correct model we can 
# obtain a plot in which "mod" is located
# inside ellipse


# Example with JSU --------------------------------------------------------

x <- rJSU(n=1000, nu=0, tau=5)
momentSK(x)

# Using the random sample
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Using the best fit from fitDist
m1 <- fitDist(x, type="realline")
m1$fits
checkMomentSK(m1)

# Using mod from gamlss function
# assuming the correct distribution
mod <- gamlss(x ~ 1, family=JSU)
mod <- refit(mod)
checkMomentSK(mod)

# Using mod from gamlss function
# assuming a wrong distribution
mod <- gamlss(x ~ 1, family=GU)
checkMomentSK(mod)

# Note: when assuming the correct model we can 
# obtain a plot in which "mod" is located
# inside ellipse
