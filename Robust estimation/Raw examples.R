
source("Robust estimation/BEor.R")
source("Robust estimation/BEorb.R")
source("Robust estimation/gamlssRobust.R")
source("Robust estimation/mywormplot.R")


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



# Creating data -----------------------------------------------------------

# Fix seed
set.seed(31416)

# sample size
number_of_obs <- 500

# proportions of mix
proportions <- sample(c(0, 1, 2), number_of_obs,
                      replace = TRUE, prob =c(0.01, 0.01, 1-0.01-0.01))

# Reference distribution (Beta(5,5))
shape_1 <- 5
shape_2 <- 5

beta_ref <- rBEo(number_of_obs, shape_1, shape_2)

# lower tail
unif_1 <- runif(number_of_obs, 0, 0.1)

# upper tail
unif_2 <- runif(number_of_obs, 0.9, 1)

#
y <- ifelse(proportions == 0, unif_1,
            ifelse(proportions == 1, unif_2, beta_ref))

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
  ylim(-0.1, 1.0) +
  xlab("")

out <- grid.arrange(p1, p2, ncol=2)
out

ggsave(plot=out, file="data_example_1.pdf", width = 6, height = 3)



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

# fitting BEo to contaminated data
# True
b1 <- gamlss(y~1, family=BEo, n.cyc=100, trace=FALSE)
p1 <- mywormplot(b1) + 
  theme_bw() + 
  ggtitle("Beta")
r1 <- cbind(fitted(b1)[1], fitted(b1, "sigma")[1])

# fit with bias correction 
b2 <- gamlss(y~1, family=BEor, n.cyc=100, trace=FALSE)
r2 <- cbind(fitted(b2)[1], fitted(b2, "sigma")[1])

#wp(b3, ylim.all=1)
p2 <- mywormplot(b2) + 
  theme_bw() + 
  ggtitle("Beta with bias correction")

# fit via robust estimation based on 
b3 <- gamlssRobust(b2, bound=2.878162, CD.bound=3.208707, trace=FALSE)
r3 <- cbind(fitted(b3)[1], fitted(b3, "sigma")[1])
#wp(b3, ylim.all=1.5)
p3 <- mywormplot(b3) + 
  theme_bw() + 
  ggtitle("Robust fitting")

out <- grid.arrange(p1, p2, p3, ncol=3)
return(out)

ggsave(plot=out, file="wp_example_1.pdf", width = 13, height = 4)


# Para crear la tabla con los parametros estimados ------------------------

true_mu    <- 5
true_sigma <- 5

Model <- c("Beta", 
           "Beta with bias correction", 
           "Robust fitting", "True")
mu <- c(fitted(b1,"mu")[1], fitted(b2,"mu")[1], 
        fitted(b3,"mu")[1], true_mu)
sigma <- c(fitted(b1,"sigma")[1], fitted(b2,"sigma")[1], 
           fitted(b3,"sigma")[1], true_sigma)

aic <- c(AIC(b1), AIC(b2), AIC(b3), NA)

tabla_resumen <- data.frame(Model, mu, sigma, aic)
tabla_resumen

library(xtable)
xtable(tabla_resumen)


# Comparando los resultados -----------------------------------------------

# Density
datos <- data.frame(y=y)

p1 <- ggplot(data = datos) + 
  geom_histogram(aes(x=y, y = ..density..),
                 bins = nclass.scott(y), 
                 colour="black", fill="white") +
  geom_function(aes(colour = "Beta"), 
                fun = dBEo, n = 10001, args = list(mu=mu[1], sigma=sigma[1])) +
  geom_function(aes(colour = "Beta with bias correction"), 
                fun = dBEo, n = 10001, args = list(mu=mu[2], sigma=sigma[2])) +
  geom_function(aes(colour = "Robust fitting"), 
                fun = dBEo, n = 10001, args = list(mu=mu[3], sigma=sigma[3])) +
  labs(x=expression(italic(y)), y=expression(italic(f(y)))) +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    legend.position = c(.03, 0.995),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
    )


# CDF

p2 <- ggplot(data = datos) + 
  stat_ecdf(data = as.data.frame(y), aes(x=y), geom="point", 
            colour=gray(.5), shape=1, size=0.9) + 
  geom_function(aes(colour = "Beta"), 
                fun = pBEo, n = 10001, args = list(mu=mu[1], sigma=sigma[1])) +
  geom_function(aes(colour = "Beta with bias correction"), 
                fun = pBEo, n = 10001, args = list(mu=mu[2], sigma=sigma[2])) +
  geom_function(aes(colour = "Robust fitting"), 
                fun = pBEo, n = 10001, args = list(mu=mu[3], sigma=sigma[3])) +
  labs(x=expression(italic(y)), y=expression(italic(F(y)))) +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    legend.position = c(.03, 0.995),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  )


# Hazard

haz_BEo <- function(x, mu, sigma) {
  dBEo(x, mu=mu, sigma=sigma) / pBEo(x, mu=mu, sigma=sigma, lower.tail=FALSE)
}

p3 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + 
  stat_function(aes(colour = "Beta"), 
                fun=haz_BEo, args=list(mu=mu[1], sigma=sigma[1])) + 
  stat_function(aes(colour = "Beta with bias correction"),
                fun=haz_BEo, args=list(mu=mu[2], sigma=sigma[2])) +
  stat_function(aes(colour = "Robust fitting"),
                fun=haz_BEo, args=list(mu=mu[3], sigma=sigma[3])) +
  xlim(0.01, 0.99) +
  labs(x=expression(italic(y)), y=expression(italic(h(y)))) +
  labs(color="Model")  +
  theme(
    legend.background = element_rect(fill = "transparent"),
    legend.position = c(.03, 0.995),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0)
  )


out <- grid.arrange(p1, p2, p3, ncol=3)

ggsave(plot=out, file="3plots_example_1.pdf", width = 13, height = 4)


