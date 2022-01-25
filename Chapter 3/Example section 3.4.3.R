library(gamlss)


# First example -----------------------------------------------------------

# generate the distribution
library(gamlss.tr)

gen.trun(par=c(0,100),family="TF", name="0to100", type="both")
## A truncated family of distributions from TF has been generated
## and saved under the names:
## dTF0to100 pTF0to100 qTF0to100 rTF0to100 TF0to100
## The type of truncation is both
## and the truncation parameter is 0 100

Y <- rTF0to100(n=500, mu=60, sigma=10, nu=5)
h1 <- histDist(Y, family=TF0to100, nbins=20, xlim=c(0, 100),
               line.col=gray(.2), line.wd=2.5, main="(a)")

# To recover the original parameters from h1
coef(h1, what='mu')
exp(coef(h1, what='sigma')) # link log
exp(coef(h1, what='nu'))    # link log


# Second example ----------------------------------------------------------

gen.trun(par=101, family="NBI", name="0to100", type="right")
## A truncated family of distributions from NBI has been generated
## and saved under the names:
## dNBI0to100 pNBI0to100 qNBI0to100 rNBI0to100 NBI0to100
## The type of truncation is right
## and the truncation parameter is 101

y <- rNBI0to100(n=500, mu=60, sigma=0.1)
h1 <- histDist(y, family=NBI0to100, nbins=100, xlim=c(0, 100),
               line.col=gray(.2),line.wd=2.5, main= "(b)")

# To recover the original parameters from h1
exp(coef(h1, what='mu'))       # link log
exp(coef(h1, what='sigma'))    # link log
