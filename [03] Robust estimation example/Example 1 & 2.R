# This script contains the R codes to replicate the examples
# of gamlssRobust( ) function with contaminated data.

# Required libraries
library(gamlss)
library(tibble)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(tidyverse)
library(CPHshape)
library(outliers)
library(EnvStats)
library(envoutliers)

# First we need to load some scripts
source("[03] Robust estimation example/BEor.R")
source("[03] Robust estimation example/BEorb.R")
source("[03] Robust estimation example/gamlssRobust.R")
source("[03] Robust estimation example/mywormplot.R")


# To simulate de datasets -------------------------------------------------

# Beta example

# Fix seed
set.seed(31416)
# sample size
number_of_obs <- 500
# proportions of mix
proportions <- sample(c(0, 1, 2), number_of_obs,
                      replace=TRUE, prob=c(0.01, 0.01, 1-0.01-0.01))
# Reference distribution (Beta(5,5))
shape_1 <- 5
shape_2 <- 5
beta_ref <- rBEo(n=number_of_obs, shape_1, shape_2)
# lower tail
unif_1 <- runif(number_of_obs, 0, 0.1)
# upper tail
unif_2 <- runif(number_of_obs, 0.9, 1)
# Simulating the y variable
y <- ifelse(proportions == 0, unif_1,
             ifelse(proportions == 1, unif_2, beta_ref))
data1 <- data.frame(y)

# ExGaussian example

# Fix seed
set.seed(31416)
# sample size
number_of_obs <- 500
# proportions of mix
proportions <- sample(c(1, 0), number_of_obs,
                      replace=TRUE, prob=c(0.9, 0.1))
# Reference ExGaussian(0.5, 0.1, 3)
true_mu    <- 0.5
true_sigma <- 0.1
true_nu    <- 3
exgauss_ref <- gamlss.dist::rexGAUS(n=number_of_obs, 
                                    mu=true_mu, 
                                    sigma=true_sigma, 
                                    nu=true_nu)
# Tail (outliers)
unift <- runif(number_of_obs, 0, 10)
# Simulating the y variable
y <- ifelse(proportions == 1, exgauss_ref, unift)
data2 <- data.frame(y)

# Histogram and boxplot ---------------------------------------------------

# Beta example

# Histogram
p1 <- ggplot(data = data1, aes(x=y)) +
  geom_histogram(aes(x=y, y=..density..), 
                 bins = nclass.scott(data1$y), 
                 colour="black", fill="white") + 
  ylab("Density")

# Boxplot
p2 <- ggplot(data1, aes(x = "", y = y)) +
  geom_boxplot(outlier.colour="tomato", outlier.shape=8,
               outlier.size=2, notch=TRUE, width = 0.4, 
               colour="black") +
  geom_jitter(colour="blue", alpha=0.2, width = 0.2) + 
  coord_flip() + 
  ylim(-0.1, 1.0) +
  xlab("")

out <- grid.arrange(p1, p2, ncol=2)

ggsave(plot=out, file="Robust estimation/data_example_1.pdf", 
       width = 6, height = 3)

# ExGaussian example

# Histogram
p1 <- ggplot(data = data2) +
  geom_histogram(aes(x=y, y=..density..), 
                 bins = nclass.scott(data2$y), 
                 colour="black", fill="white") + 
  ylab("Density")

# Boxplot
p2 <- ggplot(data2, aes(x = "", y = y)) +
  geom_boxplot(outlier.colour="tomato", outlier.shape=8,
               outlier.size=2, notch=TRUE, width = 0.4, 
               colour="black") +
  geom_jitter(colour="blue", alpha=0.2, width = 0.2) + 
  coord_flip() + 
  ylim(-0.1, 25) +
  xlab("")

out <- grid.arrange(p1, p2, ncol=2)

ggsave(plot=out, file="Robust estimation/data_example_2.pdf", 
       width = 6, height = 3)


# Fitting the models ------------------------------------------------------

# Beta example

# BEo
mod_BEo <- gamlss(y~1, family=BEo, n.cyc=100, trace=FALSE, data=data1)

p1 <- mywormplot(mod_BEo) + theme_bw() + ggtitle("Beta")
r1 <- cbind(fitted(mod_BEo)[1], 
            fitted(mod_BEo, "sigma")[1])

# fit with bias correction 
mod_BEo_correct <- gamlss(y~1, family=BEor, n.cyc=100, trace=FALSE, data=data1)

p2 <- mywormplot(mod_BEo_correct) + theme_bw() +  
  ggtitle("Beta with bias correction")
r2 <- cbind(fitted(mod_BEo_correct)[1], 
            fitted(mod_BEo_correct, "sigma")[1])

# fit via robust estimation based on 
mod_BEo_rob <- gamlssRobust(mod_BEo_correct, bound=2.878162, 
                            CD.bound=3.208707, trace=FALSE)

p3 <- mywormplot(mod_BEo_rob) + theme_bw() + 
  ggtitle("Robust fitting")
r3 <- cbind(fitted(mod_BEo_rob)[1], 
            fitted(mod_BEo_rob, "sigma")[1])

# The results

model_beta <- c("Beta", 
                "Beta with bias correction", 
                "Robust fitting", 
                "True")

mu_beta <- c(fitted(mod_BEo,"mu")[1], 
             fitted(mod_BEo_correct,"mu")[1], 
             fitted(mod_BEo_rob,"mu")[1], 
             5)

sigma_beta <- c(fitted(mod_BEo,"sigma")[1], 
                fitted(mod_BEo_correct,"sigma")[1], 
                fitted(mod_BEo_rob,"sigma")[1], 
                5)

# exGaussian example

# ExGaussian
mod_exgaus <- gamlss(y~1,family=exGAUS,n.cyc=100, trace=FALSE, data=data2)

p4 <- mywormplot(mod_exgaus)+theme_bw()+ggtitle("ExGaussian")
r4 <- cbind(fitted(mod_exgaus)[1], 
            fitted(mod_exgaus,"sigma")[1],
            fitted(mod_exgaus,"nu")[1])

# ExGaussian fit via robust estimation based  
mod_exgaus_rob <- gamlssRobust(mod_exgaus,
                               bound=2.878162, CD.bound=3.208707, trace=FALSE)

p5 <- mywormplot(mod_exgaus_rob)+theme_bw()+ggtitle("Robust ExGaussian")
r5 <- cbind(fitted(mod_exgaus_rob)[1], 
            fitted(mod_exgaus_rob,"sigma")[1],
            fitted(mod_exgaus_rob,"nu")[1]) 


# The results
model_exg <- c("ExGaussian", "Robust ExGaussian", "True")

mu_exg <- c(fitted(mod_exgaus,"mu")[1], 
            fitted(mod_exgaus_rob,"mu")[1], 
            0.5)

sigma_exg <- c(fitted(mod_exgaus,"sigma")[1], 
               fitted(mod_exgaus_rob,"sigma")[1], 
               0.1)

nu_exg <- c(fitted(mod_exgaus,"nu")[1], 
            fitted(mod_exgaus_rob,"nu")[1], 
            3)

# wormplots

out <- grid.arrange(p1, p2, p3, p4, p5, ncol=3)
return(out)

ggsave(plot=out, file="Robust estimation/wp_both_examples.pdf", 
       width = 12, height = 7)


# Comparing the results ---------------------------------------------------

# Beta example

p1_ex1 <- ggplot(data = data1) + 
  geom_histogram(aes(x=y, y = ..density..), bins=14, 
                 colour="black", fill="white") +
  geom_function(aes(colour = "Non-robust"), 
                fun = dBEo, n = 10001, 
                args = list(mu=mu_beta[1], sigma=sigma_beta[1])) +
  geom_function(aes(colour = "Robust"), 
                fun = dBEo, n = 10001, 
                args = list(mu=mu_beta[3], sigma=sigma_beta[3])) +
  geom_function(aes(colour = "True"), 
                fun = dBEo, n = 10001, 
                args = list(mu=mu_beta[4], sigma=sigma_beta[4])) +
  labs(x=expression(italic(y)), y=expression(italic(f(y))),
       title="Beta - PDF") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    legend.position = c(.03, 0.995),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_colour_manual(values = c("red", "green3", gray(.5)))

p2_ex1 <- ggplot(data = data1) + 
  stat_ecdf(aes(x=y), geom="point", 
            colour=gray(.5), shape=1, size=0.9) + 
  geom_function(aes(colour = "Non-robust"), 
                fun = pBEo, n = 10001, 
                args = list(mu=mu_beta[1], sigma=sigma_beta[1])) +
  geom_function(aes(colour = "Robust"), 
                fun = pBEo, n = 10001, 
                args = list(mu=mu_beta[3], sigma=sigma_beta[3])) +
  geom_function(aes(colour = "True"), 
                fun = pBEo, n = 10001, 
                args = list(mu=mu_beta[4], sigma=sigma_beta[4])) +
  labs(x=expression(italic(y)), y=expression(italic(F(y))),
       title="Beta - CDF") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    #legend.position = c(.03, 0.995),
    legend.position = "none",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_colour_manual(values = c("red", "green3", gray(.5)))

haz_BEo <- function(x, mu, sigma) {
  dBEo(x, mu=mu, sigma=sigma) / pBEo(x, mu=mu, sigma=sigma, lower.tail=FALSE)
}

data1$h <- haz_BEo(x=data1$y, mu=mu_beta[4], sigma=sigma_beta[4])

p3_ex1 <- ggplot(data=data1, aes(x=y, y=h)) + 
  geom_point(colour=gray(.5), shape=1, size=0.9) +
  stat_function(aes(colour = "Non-robust"), 
                fun=haz_BEo, args=list(mu=mu_beta[1], sigma=sigma_beta[1])) +
  stat_function(aes(colour = "Robust"),
              fun=haz_BEo, args=list(mu=mu_beta[3], sigma=sigma_beta[3])) +
  stat_function(aes(colour = "True"),
                fun=haz_BEo, args=list(mu=mu_beta[4], sigma=sigma_beta[4])) +
  xlim(0.01, 0.99) +
  labs(x=expression(italic(y)), y=expression(italic(h(y))),
       title="Beta - Hazard") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    #legend.position = c(.03, 0.995),
    legend.position = "none",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_colour_manual(values = c("red", "green3", gray(.5)))


# exGaussian example

p1_ex2 <- ggplot(data = data2) + 
  geom_histogram(aes(x=y, y = ..density..), bins=20,
                 colour="black", fill="white") +
  geom_function(aes(colour = "Non-robust"), 
                fun = dexGAUS, n = 10001, 
                args = list(mu=mu_exg[1], sigma=sigma_exg[1], nu=nu_exg[1])) +
  geom_function(aes(colour = "Robust"), 
                fun = dexGAUS, n = 10001, 
                args = list(mu=mu_exg[2], sigma=sigma_exg[2], nu=nu_exg[2])) +
  geom_function(aes(colour = "True"), 
                fun = dexGAUS, n = 10001, 
                args = list(mu=mu_exg[3], sigma=sigma_exg[3], nu=nu_exg[3])) +
  xlim(0.01, 20) +
  labs(x=expression(italic(y)), y=expression(italic(f(y))),
       title="exGaussian - PDF") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    legend.position = c(0.7, 1.0),
    #legend.position = "none",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_colour_manual(values = c("red", "green3", gray(.5)))


# CDF

p2_ex2 <- ggplot(data = data2) + 
  stat_ecdf(aes(x=y), geom="point", 
            colour=gray(.5), shape=1, size=0.9) + 
  geom_function(aes(colour = "Non-robust"), 
                fun = pexGAUS, n = 10001, 
                args = list(mu=mu_exg[1], sigma=sigma_exg[1], nu=nu_exg[1])) +
  geom_function(aes(colour = "Robust"), 
                fun = pexGAUS, n = 10001, 
                args = list(mu=mu_exg[2], sigma=sigma_exg[2], nu=nu_exg[2])) +
  geom_function(aes(colour = "True"), 
                fun = pexGAUS, n = 10001, 
                args = list(mu=mu_exg[3], sigma=sigma_exg[3], nu=nu_exg[3])) +
  xlim(0.01, 20) +
  labs(x=expression(italic(y)), y=expression(italic(F(y))),
       title="exGaussian - CDF") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    #legend.position = c(.8, 0.35),
    legend.position = "none",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_colour_manual(values = c("red", "green3", gray(.5)))


# Hazard

haz_exGAUS <- function(x, mu, sigma, nu) {
  dexGAUS(x, mu=mu, sigma=sigma, nu=nu) / 
    pexGAUS(x, mu=mu, sigma=sigma, nu=nu, lower.tail=FALSE)
}

data2$h <- haz_exGAUS(x=data2$y, mu=mu_exg[3], 
                       sigma=sigma_exg[3], nu=nu_exg[3])

p3_ex2 <- ggplot(data=data2, aes(x=y, y=h)) + 
  geom_point(colour=gray(.5), shape=1, size=0.9) +
  stat_function(aes(colour = "Non-robust"), 
                fun=haz_exGAUS,
                args=list(mu=mu_exg[1], sigma=sigma_exg[1], nu=nu_exg[1])) + 
  stat_function(aes(colour = "Robust"),
                fun=haz_exGAUS,
                args=list(mu=mu_exg[2], sigma=sigma_exg[2], nu=nu_exg[2])) +
  stat_function(aes(colour = "True"),
                fun=haz_exGAUS,
                args=list(mu=mu_exg[3], sigma=sigma_exg[3], nu=nu_exg[3])) +
  xlim(0.01, 1.5) +
  labs(x=expression(italic(y)), y=expression(italic(h(y))),
       title="exGaussian - Hazard") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    #legend.position = c(0.8, 0.5),
    legend.position = "red",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  ) +
  scale_colour_manual(values = c("red", "green3", gray(.5)))


# To create the joint figure

out_both_examples <- grid.arrange(p1_ex1, p2_ex1, p3_ex1, 
                                  p1_ex2, p2_ex2, p3_ex2,
                                  ncol=3)

ggsave(plot=out_both_examples, file="Robust estimation/3plots_both_examples.pdf", 
       width = 13, height = 8)
