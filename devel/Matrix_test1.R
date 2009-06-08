
library("Matrix")
library("splines")

n <- 10000
by <- 0.01
x <- runif(n)
y <- x^2 + rnorm(n, sd = 0.1)
w <- rep(1, n)
plot(x, y)

ox <- order(x)
X <- bs(x[ox], knots = seq(from = 0.1, to = 0.9, by = by),
        Boundary = c(0, 1))
class(X) <- "matrix"
dim(X)

object.size(X)

XX <- crossprod(X * w, X)
object.size(XX)
system.time(beta <- solve(XX, crossprod(X * w, y[ox])))

system.time(yhat <- X %*% beta)


sX <- as(X, "CsparseMatrix")
object.size(sX)

sXX <- crossprod(sX * w, sX)
object.size(sXX)
system.time(sbeta <- solve(sXX, crossprod(sX * w, y[ox])))

system.time(syhat <- sX %*% beta)


lines(sort(x), yhat, col = "red")
lines(sort(x), syhat + 0.1, col = "blue")
