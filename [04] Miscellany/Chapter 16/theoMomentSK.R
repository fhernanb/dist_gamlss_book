# New function theoMomentSK -----------------------------------------------
# This function can be used to obtain theoretical moments S and K
# using numerical integration inside.

# Arguments
# fam: corresponds to a gamlss family, by default is NO
# lower: minimum value of the random variable
# upper: maximum value of the random variable

theoMomentSK <- function(fam="NO", lower=-Inf, upper=Inf, ...) {
  fam <- as.gamlss.family(fam)
  fname <- fam$family[[1]]
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  xk_pdf <- function(x, k=1, ...) x^k * pdf(x, ...)
  # Moments about zero
  if (fam$type == "Discrete") {
    mu1_p <- try(add(f=xk_pdf, k=1, lower=lower, upper=upper, ...)$value, silent=TRUE)
    mu2_p <- try(add(f=xk_pdf, k=2, lower=lower, upper=upper, ...)$value, silent=TRUE)
    mu3_p <- try(add(f=xk_pdf, k=3, lower=lower, upper=upper, ...)$value, silent=TRUE)
    mu4_p <- try(add(f=xk_pdf, k=4, lower=lower, upper=upper, ...)$value, silent=TRUE)
  }
  else {
    mu1_p <- try(integrate(f=xk_pdf, k=1, lower=lower, upper=upper, ...)$value, silent=TRUE)
    mu2_p <- try(integrate(f=xk_pdf, k=2, lower=lower, upper=upper, ...)$value, silent=TRUE)
    mu3_p <- try(integrate(f=xk_pdf, k=3, lower=lower, upper=upper, ...)$value, silent=TRUE)
    mu4_p <- try(integrate(f=xk_pdf, k=4, lower=lower, upper=upper, ...)$value, silent=TRUE)
  }
  # If there are some errors in the integrals
  if (class(mu1_p) == "try-error") mu1_p <- NA
  if (class(mu2_p) == "try-error") mu2_p <- NA
  if (class(mu3_p) == "try-error") mu3_p <- NA
  if (class(mu4_p) == "try-error") mu4_p <- NA
  # Central moments
  mu2 <- mu2_p - mu1_p^2
  mu3 <- mu3_p - 3 * mu2_p * mu1_p + 2 * mu1_p^3
  mu4 <- mu4_p - 4 * mu3_p * mu1_p + 6 * mu2_p * mu1_p^2 - 3 * mu1_p^4
  # Skewness and Kurtosis
  gamma.1 <- mu3/mu2^1.5
  beta.2  <- mu4/mu2^2
  gamma.2 <- mu4/mu2^2 - 3
  # Transforming Skewness and Kurtosis
  tskew <- gamma.1/(1 + abs(gamma.1))
  tkurt <- gamma.2/(1 + abs(gamma.2))
  
  list(mom.skew = gamma.1, trans.mom.skew = tskew, 
       mom.kurt = beta.2, excess.mom.kurt = gamma.2, trans.mom.kurt = tkurt,
       expected_value = mu1_p, variance = mu2)
}


###########################################################################
#                            EXAMPLES 
###########################################################################

# Testing the function with some known distributions

library(gamlss)

theoMomentSK()
theoMomentSK(fam=NO, mu=3)
theoMomentSK(fam=NO, mu=2, sigma=2)
theoMomentSK(fam=GA, lower=0)
theoMomentSK(fam=GA, mu=2, lower=0)
theoMomentSK(fam=GA, mu=3, sigma=0.5, lower=0)
theoMomentSK(fam="TF", mu=0, sigma=1, nu=1) # Cauchy distribution does not have moments

# Using some distributions from RelDists package
library(RelDists)
theoMomentSK(fam="LIN", mu=3, lower=0)
theoMomentSK(fam="FWE", mu=2, sigma=1, lower=0)


# Comparing with momentSK using a random sample

# Gumbel distribution explained in section 18.2.1
curve(dGU(x, mu=1.5, sigma=1.7), from=0, to=20)
x <- rGU(n=1000, mu=1.5, sigma=1.7)
momentSK(x)
theoMomentSK(fam=GU, mu=1.5, sigma=1.7, lower=-Inf, upper=Inf)

# Logistic distribution explained in section 18.2.2
curve(dLO(x, mu=1.5, sigma=1.7), from=0, to=20)
x <- rLO(n=1000, mu=1.5, sigma=1.7)
momentSK(x)
theoMomentSK(fam=LO, mu=1.5, sigma=1.7, lower=-Inf, upper=Inf)

# Exponential distribution explained in section 19.2.1
curve(dEXP(x, mu=1.5), from=0, to=20)
x <- rEXP(n=1000, mu=1.5)
momentSK(x)
theoMomentSK(fam=EXP, mu=1.5, lower=0, upper=Inf)

# Inverse gaussian distribution explained in section 19.3.3
curve(dIG(x, mu=1.5, sigma=1.7), from=0, to=20)
x <- rIG(n=1000, mu=1.5, sigma=1.7)
momentSK(x)
theoMomentSK(fam=IG, mu=1.5, sigma=1.7, lower=0, upper=Inf)

# GIG distribution explained in section 19.4.4
curve(dGIG(x, mu=1, sigma=2, nu=-3), from=0, to=20)
x <- rGIG(n=1000, mu=1, sigma=2, nu=-3)
momentSK(x)
theoMomentSK(fam=GIG, mu=1, sigma=2, nu=-3, lower=0, upper=Inf)

# GB2 distribution explained in section 19.5.3
curve(dGB2(x, mu=1, sigma=2, nu=3, tau=4), from=0, to=20)
x <- rGB2(n=1000, mu=1, sigma=2, nu=3, tau=4)
momentSK(x)
theoMomentSK(fam=GB2, mu=1, sigma=2, nu=3, tau=4, lower=0, upper=Inf)


