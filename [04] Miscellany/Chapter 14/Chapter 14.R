# Example 1
library(gamlss)

sigmas <- c(0.5, 1, 2, 3.5, 6, 20)

curve(dWEI(x, mu=1, sigma=sigmas[1]), lty=1, from=0, to=3)
for (i in 2:length(sigmas)) 
  curve(dWEI(x, mu=1, sigma=sigmas[i]), add=TRUE, lty=i)
legend("topright", legend=sigmas, lty=1:length(sigmas))

p <- (1:50) / 100
theoCentileSK(fam=WEI, p=p)
