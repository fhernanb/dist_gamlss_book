
add <- function (f, lower, upper, ..., abs.tol = .Machine$double.eps^0.25) {
  
  f <- match.fun(f)
  ff <- function(x) f(x, ...)
  
  if (lower >= upper) 
    stop("invalid parameter values")
  stopifnot(length(lower) == 1, length(upper) == 1)
  
  # My auxiliar functions ------------------------------------------------------
  
  # First function
  add_minusinf_to_inf <- function(ff, ..., abs.tol) {
    x <- 0
    ans <- 0
    while (TRUE) {
      next_term <- ff(x) + ff(-x)
      ans <- ans + next_term
      if (abs(next_term) < abs.tol) break
      x <- x + 1
    }
    ans
  }
  
  # Second function
  add_lower_to_inf <- function(f, lower, ..., abs.tol) {
    x <- lower
    ans <- 0
    while (TRUE) {
      next_term <- ff(x)
      ans <- ans + next_term
      if (abs(next_term) < abs.tol) break
      x <- x + 1
    }
    ans
  }
  
  # Third function
  add_minusinf_to_upper <- function(f, upper, ..., abs.tol) {
    x <- upper
    ans <- 0
    while (TRUE) {
      next_term <- ff(x)
      ans <- ans + next_term
      if (abs(next_term) < abs.tol) break
      x <- x - 1
    }
    ans
  }
  
  # End my auxiliar functions --------------------------------------------------
  
  # Sum with finite lower and upper
  if (is.finite(lower) && is.finite(upper)) {
    wk <- sum(ff(seq(from=lower, to=upper, by=1)))
  }
  
  else {
    if (is.na(lower) || is.na(upper))
      stop("a limit is NA or NaN")
    
    if (is.finite(lower)) {
      wk <- add_lower_to_inf(ff, lower, ..., abs.tol=abs.tol)
    }
    else if (is.finite(upper)) {
      wk <- add_minusinf_to_upper(ff, upper, ..., abs.tol=abs.tol)
    }
    else {
      wk <- add_minusinf_to_inf(ff, ..., abs.tol=abs.tol)
    }
  }
  return(wk)
}

# Examples

# Poisson expected value
add(f=function(x, lambda) x*dpois(x, lambda), lower=0, upper=Inf, 
    lambda=7.5)

# Binomial expected value
add(f=function(x, size, prob) x*dbinom(x, size, prob), lower=0, upper=20, 
    size=20, prob=0.5)

# Examples with infinite series
add(f=function(x) 0.5^x, lower=0, upper=100) # Ans=2
add(f=function(x) (1/3)^(x-1), lower=1, upper=Inf) # Ans=1.5
add(f=function(x) 4/(x^2+3*x+2), lower=0, upper=Inf) # Ans=4.0
add(f=function(x) 1/(x*(log(x)^2)), lower=2, upper=Inf, abs.tol=0.000001) # Ans=2.02
add(f=function(x) 3*0.7^(x-1), lower=1, upper=Inf) # Ans=10

# Examples with wrong arguments
add(f=function(x) 0.5^x, lower=5, upper=-9)
add(f=function(x) 0.5^x, lower=5, upper=c(9, 15, 20))


