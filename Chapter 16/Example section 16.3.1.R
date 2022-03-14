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



# New function theoMomentSK -----------------------------------------------

my_theoMomentSK <- function(fam="NO", 
                            lower=-Inf, upper=Inf,
                            ...) {
  fam <- as.gamlss.family(fam)
  fname <- fam$family[[1]]
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  xk_pdf <- function(x, k=1, ...) x^k * pdf(x, ...)
  mu1_prime <- integrate(f=xk_pdf, k=1, lower=lower, upper=upper, ...)$value
  mu2_prime <- integrate(f=xk_pdf, k=2, lower=lower, upper=upper, ...)$value
  mu3_prime <- integrate(f=xk_pdf, k=3, lower=lower, upper=upper, ...)$value
  variance <- mu2_prime-(mu1_prime)^2
  rho1 <- 
  list(exp_val=mu1_prime, 
       variance=variance,
       rho1=)
}

my_theoMomentSK()
my_theoMomentSK(fam=NO, mu=3)
my_theoMomentSK(fam=NO, mu=2, sigma=2)
my_theoMomentSK(fam=GA, lower=0)
my_theoMomentSK(fam=GA, mu=2, lower=0)
my_theoMomentSK(fam=GA, mu=3, sigma=0.5, lower=0)

my_theoMomentSK(fam="TF", nu=1)


library(RelDists)
my_theoMomentSK(fam="LIN", mu=3, lower=0)
my_theoMomentSK(fam="FWE", mu=2, sigma=1, lower=0)


