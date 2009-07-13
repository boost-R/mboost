
source("bl.R")
source("helpers.R")
source("bolscw.R")
library("MASS")
library("Matrix")

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
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)
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
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)
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
