
f <- list.files(path = "../R", pattern = "R$", full = TRUE)
sapply(f, source)
library("MASS")
library("Matrix")
library("splines")

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
m2 <- fit(dpp(bols3(xn), w), y)
testfun(m1, m2)

### numeric x without intercept
m1 <- lm(y ~ xn - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols3(xn, intercept = FALSE), w), y)
testfun(m1, m2)

### factor x with intercept
m1 <- lm(y ~ xf, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols3(xf), w), y)
testfun(m1, m2)

### factor x without intercept
m1 <- lm(y ~ xf - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols3(xf, intercept = FALSE), w), y)
testfun(m1, m2)

### contrasts
m1 <- lm(y ~ xf, weights = w, contrasts = list(xf = "contr.sum"), na.action = na.exclude)
m2 <- fit(dpp(bols3(xf, contrasts.arg = list(xf = "contr.sum")), w), y)
testfun(m1, m2)

### multiple x
m1 <- lm(y ~ xn + xf, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols3(xn, xf), w), y)
testfun(m1, m2)

### interaction with binary factor
xtmp <- (z1 == "2") * xn
m1 <- lm(y ~ xtmp - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols3(xn, z = z1, intercept = FALSE), w), y)
testfun(m1, m2)

### interaction with numeric variable
m1 <- lm(y ~ z2 + z2:xn - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols3(xn, z = z2), w), y)
testfun(m1, m2)

### ridge
one <- rep(1, n)
cf1 <- coef(lm.ridge(y ~ one + xn - 1, lambda = 2))
cf2 <- coef(fit(dpp(bols3(one, xn, lambda = 2, intercept = FALSE), rep(1, n)), y))
max(abs(cf1 - cf2))
cf1 <- coef(lm.ridge(y ~ xf - 1, lambda = 2))
cf2 <- coef(fit(dpp(bols3(xf, lambda = 2, intercept = FALSE), rep(1, n)), y))
max(abs(cf1 - cf2))

### matrix (here with missing values)
cf1 <- coef(mod <- lm(y ~ xn + xf * z1 - 1, weights = w, y = TRUE, x = TRUE))
tX <- mod$x
tw <- mod$weights
ty <- mod$y
cf2 <- coef(fit(dpp(bols3(tX), weights = tw), ty))
max(abs(cf1 - cf2))

### ridge again
tX <- matrix(runif(1000), ncol = 10)
ty <- rnorm(100)
tw <- rep(1, 100)
# compute & check df
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
truedf <- sum(diag(tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)))
truedf - 2
one <- rep(1, ncol(tX))
cf1 <- coef(lm.ridge(ty ~ . - 1, data = as.data.frame(tX), lambda = la))
cf2 <- coef(fit(dpp(bols3(tX, df = 2), weights = tw), ty))
max(abs(cf1 - cf2))
# I think bols3 is better and thus right
sum((ty - tX %*% cf1)^2) + la * sum(cf1^2)
sum((ty - tX %*% cf2)^2) + la * sum(cf2^2)

# check df with weights
tw <- rpois(100, 2)
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
truedf <- sum(diag(tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)))
truedf - 2

### componentwise
cf2 <- coef(fit(dpp(bolscw(xn), weights = w), y))
cf1 <- coef(lm(y ~ xn - 1, weights = w))
max(abs(cf1 - max(cf2)))

cf2 <- coef(fit(dpp(bolscw(xn, intercept = FALSE), weights = w), y))
cf1 <- coef(lm(y ~ xn - 1, weights = w))
max(abs(cf1 - max(cf2)))

### with center
cf2 <- coef(fit(dpp(bolscw(xn, center = TRUE), weights = w), y))
cf1 <- coef(lm(y ~ 1, weights = w, subset = complete.cases(xn)))
max(abs(cf1 - max(cf2)))

cf2 <- coef(fit(dpp(bolscw(xn, intercept = FALSE, center = TRUE), weights = w), y))
tx <- xn - mean(xn, na.rm = TRUE)
cf1 <- coef(lm(y ~ tx - 1, weights = w))
max(abs(cf1 - max(cf2)))

### splines
n <- 110
x <- round(sort(runif(n, min = 0, max = 10)), 3)
f <- function(x) 1 + 0.5 * x + sin(x)
y <- f(x) + rnorm(n, sd = 0.05)
w <- rpois(n, lambda = 2)
x[sample(1:n)[1:10]] <- NA

ps <- function(d) fitted(fit(dpp(bbs3(x, df = d), w), y))
ss <- function(d) fitted(smooth.spline(x, y, w = w, df = d + 1))

sapply(1:20, function(d) c(mean((ps(d) - f(x))^2, na.rm = TRUE), 
  mean((ss(d) - f(x))^2, na.rm = TRUE),
  mean((ss(d) - ps(d))^2, na.rm = TRUE)))

max(abs(fitted(lm(y ~ x, weights = w)) - ps(0)[!is.na(x)]))

### centering
y <- y[!is.na(x)]
w <- w[!is.na(x)]
x <- x[!is.na(x)]
max(abs(fitted(fit(dpp(bbs3(x, df = 4), w), y)) - 
        fitted(lm(y ~ x)) + fitted(fit(dpp(bbs3(x, df = 1, center = TRUE), w), y))))

### varying coefficients
x1 <- runif(n, max = 2)
x2 <- sort(runif(n, max = 2 * pi))
y <- sin(x2) * x1 + rnorm(n)
w <- rep(1, n)

d <- dpp(bbs3(x2, z = x1, df = 4), w)
f <- fit(d, y)
f2 <- d$predict(list(f), newdata = data.frame(x1 = 1, x2 = x2))

max(abs(sin(x2) - f2))

set.seed(29)
x1 <- runif(n, min = -3, max = 3)
x2 <- runif(n, min = -3, max = 3)
y <- dnorm(x1) * dnorm(x2)
w <- rep(1, n)

f1 <- fit(dpp(bspatial3(x1, x2, df = 10), w), y)$fitted
f2 <- attr(mboost:::bspatial(x1, x2, df = 10), "dpp")(w)$fit(y)$fitted

X1 <- get("X", env = environment(f1))
X2 <- get("X", env = environment(f2))
max(abs(X1 - X2))

K1 <- get("K", env = environment(f1))
K2 <- get("K", env = environment(f2))
max(abs(K1 - K2))

l1 <- get("lambda", env = environment(f1))
l2 <- get("lambda", env = environment(f2))
l1
l2

df2lambda(X1, dmat = K1, df = 10, weights = w)
mboost:::df2lambda(X2, dmat = K2, df = 10, weights = w)



