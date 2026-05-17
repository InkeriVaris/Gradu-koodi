# Matern-kovarianssifunktion piirto

library(rSPDE)

# Matern-kovarianssifunktion kuvaaja eri parametrin nu arvoilla
x <- seq(from = 0, to = 10, length.out = 200)
plot(x, matern.covariance(abs(x), kappa = 1/2, nu = 1, sigma = 1),
     type = "l", ylab = "Kovarianssi", xlab = "Etäisyys"
)
lines(x, matern.covariance(abs(x), kappa = 1/2, nu = 1 / 2, sigma = 1),
      type = "l", col = "maroon3")

lines(x, matern.covariance(abs(x), kappa = 1/2, nu = 1 / 4, sigma = 1),
      type = "l", col = "steelblue2")

legend("topright",legend=c("v = 1", "v = 1/2", "v = 1/4"),
       col=c("black", "maroon3","steelblue2"), lty = (1))

# Matern-kovarianssifunktion kuvaaja eri parametrin kappa arvoilla
x <- seq(from = 0, to = 10, length.out = 200)
plot(x, matern.covariance(abs(x), kappa = 2, nu = 1, sigma = 1),
     type = "l", ylab = "Kovarianssi", xlab = "Etäisyys"
)
lines(x, matern.covariance(abs(x), kappa = 1, nu = 1, sigma = 1),
      type = "l", col = "maroon3")

lines(x, matern.covariance(abs(x), kappa = 1/2, nu = 1, sigma = 1),
      type = "l", col = "steelblue2")

legend("topright",legend=c("k = 2", "k = 1", "k = 1/2"),
       col=c("black", "maroon3","steelblue2"), lty = (1))
