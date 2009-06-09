
library("mboost")

n <- 1000
x1 <- runif(n)
x2 <- runif(n)
y <- x1^2 + x2^2 + rnorm(n, sd = 0.1)

k <- 20
Rprof("a1")
a1 <- gamboost(y ~ bspatial1(x1, x2, xknots = k, yknots = k))
Rprof(NULL)

Rprof("a2")
a2 <- gamboost(y ~ bspatial(x1, x2, xknots = k, yknots = k))
Rprof(NULL)

max(abs(coef(a1)[[1]] - coef(a2)[[1]]))

object.size(a1)
object.size(a2)

system.time(b <- bspatial(x1, x2, xknots = k, yknots = k))
system.time(bb <- attr(b, "dpp")(rep(1, length(x1))))

