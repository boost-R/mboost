
library("mboost3")
attach(asNamespace("mboost3"))
library("MASS")

### dgp
n <- 20000
xn <- round(runif(n), 3)
xn[sample(1:n)[1:(n / 100)]] <- NA
xf <- gl(4, n / 4)
xf[sample(1:n)[1:(n / 100)]] <- NA
z1 <- sample(gl(2, n / 2))
z1[sample(1:n)[1:(n / 100)]] <- NA
z2 <- round(runif(n), 3)
z2[sample(1:n)[1:(n / 100)]] <- NA
w <- rpois(n, lambda = 2)
y <- 2 * xn + rnorm(n)
y[is.na(y)] <- rnorm(sum(is.na(y)))

testfun <- function(m1, m2) {
    c(max(abs(coef(m1) - coef(m2))),
      max(abs(fitted(m1) - fitted(m2)), na.rm = TRUE))
}

### numeric x with intercept
m1 <- lm(y ~ xn, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xn), w), y)
testfun(m1, m2)

### numeric x without intercept
m1 <- lm(y ~ xn - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xn, intercept = FALSE), w), y)
testfun(m1, m2)

### factor x with intercept
m1 <- lm(y ~ xf, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xf), w), y)
testfun(m1, m2)

### factor x without intercept
m1 <- lm(y ~ xf - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xf, intercept = FALSE), w), y)
testfun(m1, m2)

### contrasts
m1 <- lm(y ~ xf, weights = w, contrasts = list(xf = "contr.sum"), na.action = na.exclude)
m2 <- fit(dpp(bols(xf, contrasts.arg = list(xf = "contr.sum")), w), y)
testfun(m1, m2)

### multiple x
m1 <- lm(y ~ xn + xf, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xn, xf), w), y)
testfun(m1, m2)

### interaction with binary factor
xtmp <- (z1 == "2") * xn
m1 <- lm(y ~ xtmp - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xn, z = z1, intercept = FALSE), w), y)
testfun(m1, m2)

### interaction with numeric variable
m1 <- lm(y ~ z2 + z2:xn - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(xn, z = z2), w), y)
testfun(m1, m2)

### ridge
one <- rep(1, n)
cf1 <- coef(lm.ridge(y ~ one + xn - 1, lambda = 2))
cf2 <- coef(fit(dpp(bols(one, xn, lambda = 2, intercept = FALSE), rep(1, n)), y))
max(abs(cf1 - cf2))
cf1 <- coef(lm.ridge(y ~ xf - 1, lambda = 2))
cf2 <- coef(fit(dpp(bols(xf, lambda = 2, intercept = FALSE), rep(1, n)), y))
max(abs(cf1 - cf2))

### matrix (here with missing values)
cf1 <- coef(mod <- lm(y ~ xn + xf * z1 - 1, weights = w, y = TRUE, x = TRUE))
tX <- mod$x
tw <- mod$weights
ty <- mod$y
cf2 <- coef(fit(dpp(bols(tX), weights = tw), ty))
max(abs(cf1 - cf2))

### ridge again with matrix interface
tX <- matrix(runif(1000), ncol = 10)
ty <- rnorm(100)
tw <- rep(1, 100)
# compute & check df
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
truedf <- sum(diag(tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)))
truedf - 2
one <- rep(1, ncol(tX))
cf1 <- coef(lm.ridge(ty ~ . - 1, data = as.data.frame(tX), lambda = la))
cf2 <- coef(fit(dpp(bols(tX, df = 2), weights = tw), ty))
max(abs(cf1 - cf2))
# I think bols is better and thus right
sum((ty - tX %*% cf1)^2) + la * sum(cf1^2)
sum((ty - tX %*% cf2)^2) + la * sum(cf2^2)

# check df with weights
tw <- rpois(100, 2)
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
truedf <- sum(diag(tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)))
truedf - 2

### componentwise
cf2 <- coef(fit(dpp(bolscw(cbind(1, xn)), weights = w), y))
cf1 <- coef(lm(y ~ xn - 1, weights = w))
max(abs(cf1 - max(cf2)))

cf2 <- coef(fit(dpp(bolscw(matrix(xn, nc = 1)), weights = w), y))
cf1 <- coef(lm(y ~ xn - 1, weights = w))
max(abs(cf1 - max(cf2)))

### componentwise with matrix
n <- 200
m <- 10000
x <- rnorm(n * m)
x[abs(x) < 2] <- 0
X <- Matrix(data = x, ncol = m, nrow = n)
beta <- rpois(ncol(X), lambda = 1)
y <- X %*% beta + rnorm(nrow(X))
w <- rep(1, nrow(X)) ###rpois(nrow(X), lambda = 1)
f1 <- dpp(bolscw(X), weights = w)$fit
system.time(for (i in 1:100) f <- f1(y))

### splines
n <- 110
x <- round(sort(runif(n, min = 0, max = 10)), 3)
f <- function(x) 1 + 0.5 * x + sin(x)
y <- f(x) + rnorm(n, sd = 0.05)
w <- rpois(n, lambda = 2)
x[sample(1:n)[1:10]] <- NA

ps <- function(d) fitted(fit(dpp(bbs(x, df = d), w), y))
ss <- function(d) fitted(smooth.spline(x, y, w = w, df = d + 1))

sapply(1:20, function(d) c(mean((ps(d) - f(x))^2, na.rm = TRUE), 
  mean((ss(d) - f(x))^2, na.rm = TRUE),
  mean((ss(d) - ps(d))^2, na.rm = TRUE)))

max(abs(fitted(lm(y ~ x, weights = w)) - ps(0)[!is.na(x)]))

### varying coefficients
x1 <- runif(n, max = 2)
x2 <- sort(runif(n, max = 2 * pi))
y <- sin(x2) * x1 + rnorm(n)
w <- rep(1, n)

d <- dpp(bbs(x2, z = x1, df = 4), w)
f <- fit(d, y)
f2 <- d$predict(list(f), newdata = data.frame(x1 = 1, x2 = x2))

max(abs(sin(x2) - f2))

### bols and bbs; matrix interfaces
n <- 10000
x <- runif(n, min = 0, max = 2*pi)
y <- sin(x) + rnorm(n, sd = 0.1)
w <- rpois(n, lambda = 1)
x[sample(1:n)[1:(n / 100)]] <- NA
h <- hyper_bbs(data.frame(x = x), vary = "")
X <- X_bbs(data.frame(x = x), vary = "", h)$X
f1 <- fit(dpp(bbs(x, df = ncol(X)), w), y)
f2 <- fit(dpp(bols(X, df = ncol(X)), w), y)
max(abs(coef(f1) - coef(f2)))

all.equal(get_index(data.frame(x, x)), get_index(X))
all.equal(get_index(data.frame(x)), get_index(X))

### spatial
set.seed(29)
x1 <- runif(n, min = -3, max = 3)
x2 <- runif(n, min = -3, max = 3)
y <- dnorm(x1) * dnorm(x2)
w <- rep(1, n)

f1 <- fit(dpp(bspatial(x1, x2, df = 12), w), y)$fitted
f2 <- attr(mboost:::bspatial(x1, x2, df = 12), "dpp")(w)$fit(y)$fitted

X1 <- get("X", env = environment(f1))
X2 <- get("X", env = environment(f2))
max(abs(X1 - X2))

K1 <- get("K", env = environment(f1))
K2 <- get("K", env = environment(f2))
max(abs(K1 - K2))

l1 <- get("lambda", env = environment(f1))
l2 <- get("lambda", env = environment(f2))
l1 - l2
### because of near-zero eigenvalues.

x <- seq(from = 0, to = 2 * pi, length = 100)
y <- sin(x) + rnorm(length(x), sd = 0.1)
m1 <- mboost(y ~ mboost3:::bbs(x, df = 1, center = TRUE) + mboost3:::bols(x))
m2 <- mboost:::gamboost(y ~ mboost:::bbs(x, df = 1, center = TRUE) + mboost:::bols(x))
stopifnot(max(abs(fitted(m1) - fitted(m2))) < sqrt(.Machine$double.eps))
