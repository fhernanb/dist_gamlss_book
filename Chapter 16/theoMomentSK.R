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
  mu1_prime <- integrate(f=xk_pdf, k=1, lower=lower, upper=upper, ...)$value
  mu2_prime <- integrate(f=xk_pdf, k=2, lower=lower, upper=upper, ...)$value
  mu3_prime <- integrate(f=xk_pdf, k=3, lower=lower, upper=upper, ...)$value
  mu4_prime <- integrate(f=xk_pdf, k=4, lower=lower, upper=upper, ...)$value
  mu2 <- mu2_prime - mu1_prime^2
  mu3 <- mu3_prime - 3 * mu2_prime * mu1_prime + 2 * mu1_prime^3
  mu4 <- mu4_prime - 4 * mu3_prime * mu1_prime + 6 * mu2_prime * mu1_prime^2 - 3 * mu1_prime^4
  gamma.1 <- mu3/mu2^1.5
  beta.2  <- mu4/mu2^2
  gamma.2 <- mu4/mu2^2 - 3
  tskew <- gamma.1/(1 + abs(gamma.1))
  tkurt <- gamma.2/(1 + abs(gamma.2))
  
  list(mom.skew = gamma.1, trans.mom.skew = tskew, 
       mom.kurt = beta.2, excess.mom.kurt = gamma.2, trans.mom.kurt = tkurt,
       expected_value = mu1_prime, variance = mu2)
}


# Testing the function with some known distributions

library(gamlss)
theoMomentSK()
theoMomentSK(fam=NO, mu=3)
theoMomentSK(fam=NO, mu=2, sigma=2)
theoMomentSK(fam=GA, lower=0)
theoMomentSK(fam=GA, mu=2, lower=0)
theoMomentSK(fam=GA, mu=3, sigma=0.5, lower=0)

theoMomentSK(fam="TF", nu=1)


library(RelDists)
theoMomentSK(fam="LIN", mu=3, lower=0)
theoMomentSK(fam="FWE", mu=2, sigma=1, lower=0)


# Comparing with momentSK using a random sample

# Gumbel distribution explained in section 18.2.1
x <- rGU(n=1000000, mu=1.5, sigma=1.7)
momentSK(x)
theoMomentSK(fam=GU, mu=1.5, sigma=1.7, lower=-Inf, upper=Inf)

# Logistic distribution explained in section 18.2.2
x <- rLO(n=1000000, mu=1.5, sigma=1.7)
momentSK(x)
theoMomentSK(fam=LO, mu=1.5, sigma=1.7, lower=-Inf, upper=Inf)

# Exponential distribution explained in section 19.2.1
x <- rEXP(n=1000000, mu=1.5)
momentSK(x)
theoMomentSK(fam=EXP, mu=1.5, lower=0, upper=Inf)

# Inverse gaussian distribution explained in section 19.3.3
x <- rIG(n=1000, mu=1.5, sigma=1.7)
momentSK(x)
theoMomentSK(fam=IG, mu=1.5, sigma=1.7, lower=0, upper=Inf)


