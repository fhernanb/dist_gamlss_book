
library(gamlss)

dax  <- EuStockMarkets[,"DAX"]
Rdax <- diff(log(dax))

m1 <- fitDist(Rdax)
m1$fits

momentSK(Rdax)

checkMomentSK(Rdax, col.b=gray(.6), pch.b=4); title("(a)")
checkMomentSK(m1,   col.b=gray(.6), pch.b=4); title("(b)")


# Using centile skewness and kurtosis -------------------------------------

checkCentileSK(Rdax, col.b=gray(.6), pch.b=4)
title("(a)")

checkCentileSK(Rdax, type="tail", col.b=gray(.6), pch.b=4)
title("(b)")

checkCentileSK(m1, col.b=gray(.6), pch.b=4)
title("(c)")

checkCentileSK(m1, type="tail", col.b=gray(.6), pch.b=4)
title("(d)")


# Using random sample -----------------------------------------------------

# Ejemplo NO
x <- rNO(n=1000, mu=170, sigma=10)
momentSK(x)
checkMomentSK(x, col.b=gray(.6), pch.b=4)
# Con el modelo
m1 <- fitDist(x)
m1$fits
checkMomentSK(m1)

mod <- gamlss(x ~ 1, family=NO)
checkMomentSK(mod)

mod <- gamlss(x ~ 1, family=WEI)
checkMomentSK(mod)

# Ejemplo GA
x <- rGA(n=1000, mu=170, sigma=10)
momentSK(x)
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Con el modelo
m1 <- fitDist(x)
m1$fits
checkMomentSK(m1)

# Ejemplo PO
x <- rPO(n=10000, mu=170)
momentSK(x)
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Con el modelo
m1 <- fitDist(x, type="counts")
m1$fits
checkMomentSK(m1)


# Usando simulacion de las distribuciones de checkmomentsk ----------------
x <- rJSU(n=1000, nu=0, tau=5)

plot(density(x))

momentSK(x)
checkMomentSK(x, col.b=gray(.6), pch.b=4)

# Con el modelo
m1 <- fitDist(x, type="realAll")
m1$fits
checkMomentSK(m1)


# Explorando todas las funciones ------------------------------------------
x <- rGA(n=1000, mu=170, sigma=10)

momentSK(x,  weights=NULL)
centileSK(x, cent = c(1, 25), weights=NULL)
centileSkew(x, cent = 1, weights=NULL)
centileKurt(x, cent = 1, weights=NULL)


my_theoMomentSK <- function(fam="NO", ...) {
  fam <- as.gamlss.family(fam)
  fname <- fam$family[[1]]
  dfun <- paste("d", fname, sep = "")
  pdf <- eval(parse(text = dfun))
  #curve(pdf(x, ...), from=0.01, to=10)
  #integrate(f=pdf(x, ...), lower=0, upper=2)
  #class(pdf)
  list(...)
}

my_theoMomentSK()
my_theoMomentSK(fam=GA)
my_theoMomentSK(fam="PO")

library(RelDists)
my_theoMomentSK(fam="LIN", mu=3)
my_theoMomentSK(fam="FWE", mu=2, sigma=1)


