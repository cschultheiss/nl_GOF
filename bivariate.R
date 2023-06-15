require(mgcv)
require(FOCI)
pot <- function(x, b) sign(x)*abs(x)^b
va <- function(b) 2^(b)*gamma(b + 0.5)/sqrt(pi)


fx1 <- function(x) dnorm(x)
f1 <- function(x) sin(2 * x)
b <- 1
scale.a <- sqrt(6)
fx2a <- function(x) dexp(abs(x)) / 2
fx2b <- function(x) dnorm(x, sd = 1)
fx2 <- Vectorize(function(x) integrate(function(z) fx2a(x-z)*fx2b(z),-Inf,Inf)$value)
fx2 <- function(x) dnorm(x, sd = sqrt(0.5))

fx12 <- function(x1, x2) fx1(x1) * fx2(x2 -f1(x1))

n <- 1e5
x1 <- rnorm(n)
epsa <- rnorm(n, sd = sqrt(0.5))
x2 <- f1(x1) + epsa
par(mfrow = c(1,1))
plot(fx2, xlim = range(epsa))
lines(density(epsa), col = 2)

par(mfrow = c(3,3))
for (x2o in quantile(x2, (1:9)/ 10)){
  plot(function(x1) fx12(x1, x2o), xlim = c(2 * min(x1), 2 * max(x1)))
}

for (q in (1:9)/10){
  plot(density(x1[(x2 > quantile(x2, q - 0.05)) & (x2 < quantile(x2, q + 0.05))]))
}

x.grid <- seq(quantile(x2, 0.001), quantile(x2, 0.999), length.out = 1000)

norm <- Vectorize(function(x2) integrate(function(x1) fx12(x1, x2), 2 * min(x1), 2 * max(x1))$value)
m1 <- Vectorize(function(x2) integrate(function(x1) x1 * fx12(x1, x2), 2 * min(x1), 2 * max(x1))$value)
# m1 <- function(x2) m1(x2)/norm(x2)
m2 <- Vectorize(function(x2) integrate(function(x1) x1^2 * fx12(x1, x2), 2 * min(x1), 2 * max(x1))$value)
# m2 <- function(x2) m2(x2)/norm(x2)
m3 <- Vectorize(function(x2) integrate(function(x1) x1^3 * fx12(x1, x2), 2 * min(x1), 2 * max(x1))$value)
# m3 <- function(x2) m3(x2)/norm(x2)

norm.grid <- norm(x.grid)
ss1 <- smooth.spline(x.grid, m1(x.grid)/norm.grid, cv = TRUE)
ss2 <- smooth.spline(x.grid, m2(x.grid)/norm.grid, cv = TRUE)
ss3 <- smooth.spline(x.grid, m3(x.grid)/norm.grid, cv = TRUE)

m1x <- predict(ss1, x2)$y
m2x <- predict(ss2, x2)$y
m3x <- predict(ss3, x2)$y
ga1 <- gam(x1 ~ s(x2))
ga2 <- gam(x1^2 ~ s(x2))
ga3 <- gam(x1^3 ~ s(x2))

subind <- which((x2 > quantile(x2, 0.01)) & (x2 < quantile(x2, 0.99)))
par(mfrow = c(1,1))
plot(x2[subind], (m2x- m1x^2)[subind])
abline(v = quantile(x2, c(0.001, 0.999)))
points(x2, ga2$fitted.values - ga1$fitted.values^2, col = 2)

plot(x2[subind], ((m3x - 3 * m1x * (m2x - m1x^2) - m1x^3)/(m2x - m1x^2)^1.5)[subind])

eps0u <- (x1- m1x)
eps0 <- (x1- m1x)/sqrt(m2x - m1x^2)
epsu <- (x1 - ga1$fitted.values)
eps <- (x1 - ga1$fitted.values)/sqrt(ga2$fitted.values - ga1$fitted.values^2)
cor(eps0, eps, use ="complete.obs", method = "spearman")

plot(eps0[subind], eps[subind])
codec(eps0u, x2)
codec(epsu, x2)
codec(eps0, x2)
codec(eps, x2)
codec(abs(eps0) + sign(eps0), x2)
codec(abs(eps) + sign(eps), x2)
plot(x2, eps0)
