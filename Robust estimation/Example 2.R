
source("Robust estimation/gamlssRobust.R")
source("Robust estimation/mywormplot.R")


library(gamlss)
library(tibble)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(tidyverse)
library(RelDists)

library(CPHshape)
library(outliers)
library(EnvStats)
library(envoutliers)



# Creating data -----------------------------------------------------------

# Fix seed
set.seed(31416)

# sample size
number_of_obs <- 500

# proportions of mix
proportions <- sample(c(1, 0), number_of_obs,
                      replace = TRUE, prob =c(0.9, 0.1))

# Reference ExGaussian(.5,0.1,3)
true_mu    <- 0.5
true_sigma <- 0.1
true_nu    <- 3

exgauss_ref <- gamlss.dist::rexGAUS(number_of_obs, 
                                    mu=true_mu, 
                                    sigma=true_sigma, 
                                    nu=true_nu)

# Tail (outliers)
unift <- runif(number_of_obs, 0, 10)

# response
y <- ifelse(proportions == 1, exgauss_ref, unift)

datos <- data.frame(y)


# Creando el histograma y boxplot -----------------------------------------

# Histograma
p1 <- ggplot(data = datos) +
  geom_histogram(aes(x=y, y=..density..), 
                 bins = nclass.scott(datos$y), 
                 colour="black", fill="white") + 
  ylab("Density")

# Boxplot
p2 <- ggplot(datos, aes(x = "", y = y)) +
  geom_boxplot(outlier.colour="tomato", outlier.shape=8,
               outlier.size=2, notch=TRUE, width = 0.4, 
               colour="black") +
  geom_jitter(colour="blue", alpha=0.2, width = 0.2) + 
  coord_flip() + 
  ylim(-0.1, 25) +
  xlab("")

out <- grid.arrange(p1, p2, ncol=2)
out

ggsave(plot=out, file="data_example_2.pdf", width = 6, height = 3)



# Identificando outliers --------------------------------------------------

# Detection of outliers - IQR criterion
value <- boxplot.stats(y)$out
index <- which(y %in% c(value))
outliers_1 <- tibble("index"=index, "value"=value)
outliers_1

# Detection of outliers - Hampel

# Bounds of Hampel criterion
lower_bound <- median(y) - 3 * mad(y, constant = 1)
upper_bound <- median(y) + 3 * mad(y, constant = 1)
# identification
index <- which(y < lower_bound | y > upper_bound)
value <- y[index]
outliers_2 <- tibble("index"=index, "value"=value)
outliers_2

# Detection of outliers - Rosner test-k
out_3 <- rosnerTest(y, k=3, warn=FALSE)$all.stats
out_3 <- out_3[order(out_3[,5]),]
as_tibble(out_3)

out_14 <- rosnerTest(y, k=14, warn=FALSE)$all.stats
out_14 <- out_14[order(out_14[,5]),]
as_tibble(out_14)


# Para crear la tabla de outliers

temp <- outliers_2
temp[1:14, 1:2] <- NA
temp[1:3, 1:2] <- outliers_1
temp

los_out <- cbind(temp, outliers_2)

xtable(los_out)

# Ajustando los modelos ---------------------------------------------------

# fitting Exgaussian to contaminated data

# True model
true_mod <- gamlss(y~1,family=exGAUS,n.cyc=100, trace=FALSE,
                   mu.fix = TRUE, mu.start=true_mu,
                   sigma.fix = TRUE, sigma.start=true_sigma,
                   nu.fix=TRUE, nu.start=true_nu)

p.0 <- mywormplot(true_mod)+theme_bw()+ggtitle("True model")

# ExGaussian
mod_exgaus <- gamlss(y~1,family=exGAUS,n.cyc=100, trace=FALSE)
p.1 <- mywormplot(mod_exgaus)+theme_bw()+ggtitle("ExGaussian")
r.1 <- cbind(fitted(mod_exgaus)[1], fitted(mod_exgaus,"sigma")[1],fitted(mod_exgaus,"nu")[1])

# ExGaussian fit via robust estimation based  
mod_exgaus_rob <- gamlssRobust(mod_exgaus,
                               bound=2.878162, CD.bound=3.208707,trace=FALSE)
r.2 <- cbind(fitted(mod_exgaus_rob)[1], fitted(mod_exgaus_rob,"sigma")[1],fitted(mod_exgaus_rob,"nu")[1]) 
p.2 <- mywormplot(mod_exgaus_rob)+theme_bw()+ggtitle("Robust ExGaussian")

# Lindley - Based on RelDists package
mod_lin <- gamlss(y~1,family=LIN,n.cyc=100, trace=FALSE)
p.3 <- mywormplot(mod_lin)+theme_bw()+ggtitle("Lindley")
r.3 <- cbind(fitted(mod_lin)[1])

# Gamma
mod_gam <- gamlss(y~1,family=GA,n.cyc=100, trace=FALSE)
p.4 <- mywormplot(mod_gam)+theme_bw()+ggtitle("Gamma")
r.4 <- cbind(fitted(mod_gam)[1], fitted(mod_gam,"sigma")[1])

# IG - Inverse Gaussian
mod_invgaus <- gamlss(y~1,family=IG,n.cyc=100, trace=FALSE)
p.5 <- mywormplot(mod_invgaus)+theme_bw()+ggtitle("Inverse Gaussian")
r.5 <- cbind(fitted(mod_invgaus)[1], fitted(mod_invgaus,"sigma")[1])

# Weibull
mod_wei <- gamlss(y~1,family=IG,n.cyc=100, trace=FALSE)
p.6 <- mywormplot(mod_wei)+theme_bw()+ggtitle("Weibull")
r.6 <- cbind(fitted(mod_wei)[1], fitted(mod_wei,"sigma")[1])

out <- grid.arrange(p.1, p.2, p.3, p.4, p.5, p.6, ncol=3)
out

ggsave(plot=out, file="wp_example_2.pdf", width = 13, height = 8)

ggsave(plot=grid.arrange(p.0, p.1, p.2, ncol=3), 
       file="wp_example_2_true_usual_robust.pdf", 
       width = 13, height = 4)



# Para crear la tabla con los parametros estimados ------------------------

true_mu    <- 0.5
true_sigma <- 0.1
true_nu    <- 3

Model <- c("ExGaussian", "Robust ExGaussian", 
           "Lindley", "Gamma", "Inverse Gaussian", "Weibull",
           "True")
mu <- c(fitted(mod_exgaus,"mu")[1], fitted(mod_exgaus_rob,"mu")[1], 
        fitted(mod_lin,"mu")[1], fitted(mod_gam,"mu")[1],
        fitted(mod_invgaus,"mu")[1], fitted(mod_wei,"mu")[1],
        true_mu)
sigma <- c(fitted(mod_exgaus,"sigma")[1], fitted(mod_exgaus_rob,"sigma")[1], 
           NA, fitted(mod_gam,"sigma")[1],
           fitted(mod_invgaus,"sigma")[1], fitted(mod_wei,"sigma")[1],
           true_sigma)
nu <- c(fitted(mod_exgaus,"nu")[1], fitted(mod_exgaus_rob,"nu")[1], 
        NA, NA, NA, NA,
        true_nu)

aic <- c(AIC(mod_exgaus), AIC(mod_exgaus_rob), AIC(mod_lin), AIC(mod_gam), AIC(mod_invgaus), AIC(mod_wei), NA)

tabla_resumen <- data.frame(Model, mu, sigma, nu, aic)
tabla_resumen

library(xtable)
xtable(tabla_resumen)


# Comparando los resultados -----------------------------------------------

# Density
datos <- data.frame(y=y)

p1_ex2 <- ggplot(data = datos) + 
  geom_histogram(aes(x=y, y = ..density..),
                 bins = nclass.scott(y), 
                 colour="black", fill="white") +
  geom_function(aes(colour = "Non-robust"), 
                fun = dexGAUS, n = 10001, 
                args = list(mu=mu[1], sigma=sigma[1], nu=nu[1])) +
  geom_function(aes(colour = "Robust"), 
                fun = dexGAUS, n = 10001, 
                args = list(mu=mu[2], sigma=sigma[2], nu=nu[2])) +
  geom_function(aes(colour = "True"), 
                fun = dexGAUS, n = 10001, 
                args = list(mu=mu[7], sigma=sigma[7], nu=nu[7])) +
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
    )


# CDF

p2_ex2 <- ggplot(data = datos) + 
  stat_ecdf(data = as.data.frame(y), aes(x=y), geom="point", 
            colour=gray(.5), shape=1, size=0.9) + 
  geom_function(aes(colour = "Non-robust"), 
                fun = pexGAUS, n = 10001, 
                args = list(mu=mu[1], sigma=sigma[1], nu=nu[1])) +
  geom_function(aes(colour = "Robust"), 
                fun = pexGAUS, n = 10001, 
                args = list(mu=mu[2], sigma=sigma[2], nu=nu[2])) +
  geom_function(aes(colour = "True"), 
                fun = pexGAUS, n = 10001, 
                args = list(mu=mu[7], sigma=sigma[7], nu=nu[7])) +
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
  )


# Hazard

haz_exGAUS <- function(x, mu, sigma, nu) {
  dexGAUS(x, mu=mu, sigma=sigma, nu=nu) / pexGAUS(x, mu=mu, sigma=sigma, nu=nu, lower.tail=FALSE)
}

p3_ex2 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + 
  stat_function(aes(colour = "Non-robust"), 
                fun=haz_exGAUS,
                args=list(mu=mu[1], sigma=sigma[1], nu=nu[1])) + 
  stat_function(aes(colour = "Robust"),
                fun=haz_exGAUS,
                args=list(mu=mu[2], sigma=sigma[2], nu=nu[2])) +
  stat_function(aes(colour = "True"),
                fun=haz_exGAUS,
                args=list(mu=mu[7], sigma=sigma[7], nu=nu[7])) +
  xlim(0.01, 1.5) +
  labs(x=expression(italic(y)), y=expression(italic(h(y))),
       title="exGaussian - Hazard") +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    #legend.position = c(0.8, 0.5),
    legend.position = "none",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  )


out_ex2 <- grid.arrange(p1_ex2, p2_ex2, p3_ex2, ncol=3)

ggsave(plot=out_ex2, file="3plots_example_2.pdf", width = 13, height = 4)


# Para crear la figura conjunta -------------------------------------------

out_both_examples <- grid.arrange(p1_ex1, p2_ex1, p3_ex1, 
                        p1_ex2, p2_ex2, p3_ex2,
                        ncol=3)

ggsave(plot=out_both_examples, file="3plots_both_examples.pdf", 
       width = 13, height = 8)


