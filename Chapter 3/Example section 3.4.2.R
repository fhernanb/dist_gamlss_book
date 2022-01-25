library(gamlss)

# Problem: in this example it is possible to obtain y=0 and/or y=1. We
# explain how to solve this issue.

gen.Family("TF", type="logit")
## A logit family of distributions from TF has been generated
## and saved under the names:
## dlogitTF plogitTF qlogitTF rlogitTF logitTF

# Generate 200 observations with mu=0, sigma=1, nu=1.5
set.seed(1234567)
Y <- rlogitTF(200, mu=0, sigma=1, nu=1.5)

# To obtain a summary of Y
min(Y)
max(Y)

# If we use the original Y we will find errors because the logitTF distribution
# does not accept 0 or 1 values.

# To avoid the 0 and 1 values not accepted in the logitTF distribution
# we can use this code
Y[Y == 0] <- 0.0001
Y[Y == 1] <- 0.9999

# fit the distribution
h1 <- histDist(Y, family=logitTF, nbins=20, ylim=c(0,2), xlim=c(0,1),
               line.col=1, nline.wd=2.5, main="(b)")

# We can fit a model to recover the original parameters
mod <- gamlss(Y ~ 1, family='logitTF')
coef(mod, what='mu')
exp(coef(mod, what='sigma')) # link log
exp(coef(mod, what='nu'))    # link log
