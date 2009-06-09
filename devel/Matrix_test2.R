
library("mboost")

n <- 20000
x <- runif(n)
y <- x^2 + rnorm(n, sd = 0.1)

k <- 40
Rprof("a1")
a1 <- gamboost(y ~ bbs1(x, knots = k))
Rprof(NULL)

Rprof("a2")
a2 <- gamboost(y ~ bbs(x, knots = k))
Rprof(NULL)

max(abs(coef(a1)[[1]] - coef(a2)[[1]]))

object.size(a1)
object.size(a2)

