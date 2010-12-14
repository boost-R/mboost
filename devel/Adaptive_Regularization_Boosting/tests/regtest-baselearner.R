
library("mboost")
attach(asNamespace("mboost"))
library("MASS")
library("Matrix")

set.seed(290875)

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
    ret <- c(max(abs(coef(m1) - coef(m2))),
      max(abs(fitted(m1) - fitted(m2)), na.rm = TRUE))
    if (any(ret > sqrt(.Machine$double.eps)))
        return(ret)
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
m2 <- fit(dpp(bols(xn, by = z1, intercept = FALSE), w), y)
testfun(m1, m2)

### interaction with numeric variable
m1 <- lm(y ~ z2:xn - 1, weights = w, na.action = na.exclude)
m2 <- fit(dpp(bols(z2, by = xn, intercept = FALSE), w), y)
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
stopifnot(max(abs(cf1 - cf2)) < sqrt(.Machine$double.eps))

### ridge again with matrix interface
tX <- matrix(runif(1000), ncol = 10)
ty <- rnorm(100)
tw <- rep(1, 100)

### compute & check df
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
truedf <- sum(diag(tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)))
stopifnot(abs(truedf - 2) < sqrt(.Machine$double.eps))

one <- rep(1, ncol(tX))
cf1 <- coef(lm.ridge(ty ~ . - 1, data = as.data.frame(tX), lambda = la))
cf2 <- coef(fit(dpp(bols(tX, df = 2), weights = tw), ty))
max(abs(cf1 - cf2))
# I think bols is better and thus right
sum((ty - tX %*% cf1)^2) + la * sum(cf1^2)
sum((ty - tX %*% cf2)^2) + la * sum(cf2^2)

### now with other df-definition:
op <- options(mboost_dftraceS = FALSE)
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
H <- tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)
truedf <- sum(diag(2*H - tcrossprod(H,H)))
stopifnot(abs(truedf - 2) < sqrt(.Machine$double.eps))
options(op)

# check df with weights
tw <- rpois(100, 2)
la <- df2lambda(tX, df = 2, dmat = diag(ncol(tX)), weights = tw)["lambda"]
truedf <- sum(diag(tX %*% solve(crossprod(tX * tw, tX) + la * diag(ncol(tX))) %*% t(tX * tw)))
stopifnot(abs(truedf - 2) < sqrt(.Machine$double.eps))

### check df2lambda for P-splines (Bug spotted by B. Hofner)
set.seed(1907)
x <- runif(100, min = -1, max = 3)
## extract lambda from base-learner
lambda <- bbs(x, df = 4)$dpp(rep(1, length(x)))$df()["lambda"]
X <- get("X", envir = environment(bbs(x, df = 4)$dpp))
K <- get("K", envir = environment(bbs(x, df = 4)$dpp))
truedf <- sum(diag(X %*%  solve(crossprod(X,X) + lambda * K) %*% t(X)))
stopifnot(abs(truedf - 4) < sqrt(.Machine$double.eps))

### check accuracy of df2lambda
data("bodyfat", package="mboost")
diff_df <- matrix(NA, nrow=8, ncol=ncol(bodyfat))
rownames(diff_df) <- paste("df", 3:10)
colnames(diff_df) <- names(bodyfat)
for (i in 3:10){
    for (j in 1:ncol(bodyfat)){
        lambda <- bbs(bodyfat[[j]], df = i)$dpp(rep(1, nrow(bodyfat)))$df()["lambda"]
        diff_df[i-2,j] <- bbs(bodyfat[[j]], lambda = lambda)$dpp(rep(1, nrow(bodyfat)))$df()["df"] - i
    }
}
stopifnot(all(diff_df < sqrt(.Machine$double.eps)))

### componentwise
cf2 <- coef(fit(dpp(bolscw(cbind(1, xn)), weights = w), y))
cf1 <- coef(lm(y ~ xn - 1, weights = w))
stopifnot(max(abs(cf1 - max(cf2))) < sqrt(.Machine$double.eps))

cf2 <- coef(fit(dpp(bolscw(matrix(xn, nc = 1)), weights = w), y))
cf1 <- coef(lm(y ~ xn - 1, weights = w))
stopifnot(max(abs(cf1 - max(cf2))) < sqrt(.Machine$double.eps))

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
f1(y)$model

### varying coefficients
x1 <- runif(n, max = 2)
x2 <- sort(runif(n, max = 2 * pi))
y <- sin(x2) * x1 + rnorm(n)
w <- rep(1, n)

d <- dpp(bbs(x2, by = x1, df = 4), w)
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
stopifnot(max(abs(coef(f1) - coef(f2))) < sqrt(.Machine$double.eps))

stopifnot(all.equal(get_index(data.frame(x, x)), get_index(X)))
stopifnot(all.equal(get_index(data.frame(x)), get_index(X)))


### combinations and tensor products of base-learners

set.seed(29)
n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
f <- gl(4, 25)
y <- rnorm(n)
ndf <- data.frame(x1 = x1[1:10], x2 = x2[1:10], f = f[1:10])

### spatial
m1 <- gamboost(y ~ bbs(x1) %X% bbs(x2))
m2 <- gamboost(y ~ bspatial(x1, x2, df = 16))
stopifnot(max(abs(predict(m1) - predict(m2))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(m1, newdata = ndf) - predict(m2, newdata = ndf))) < sqrt(.Machine$double.eps))

### spatio-temporal
m1 <- gamboost(y ~ bbs(x1, knots = 6) %X% bbs(x2, knots = 6) %X% bbs(x3, knots = 6))
m2 <- gamboost(y ~ (bbs(x1, knots = 6) + bbs(x2, knots = 6)) %X% bbs(x3, knots = 6))

### varying numeric
m1 <- gamboost(y ~ bbs(x1) %X% bols(x2, intercept = FALSE, lambda = 0))
m2 <- gamboost(y ~ bbs(x1, by = x2, df = 4))
stopifnot(max(abs(predict(m1) - predict(m2))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(m1, newdata = ndf) - predict(m2, newdata = ndf))) < sqrt(.Machine$double.eps))

### varying factor
m1 <- gamboost(y ~ bbs(x1) %X% bols(f, intercept = FALSE, df = 5))
coef(m1)
predict(m1, newdata = ndf)

### cbind
m1 <- gamboost(y ~ bols(x1, intercept = FALSE, df = 1) %+%
                   bols(x2, intercept = FALSE, df = 1))
m2 <- gamboost(y ~ bols(x1, x2, intercept = FALSE, df = 2))
stopifnot(max(abs(predict(m1) - predict(m2))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(m1, newdata = ndf) - predict(m2, newdata = ndf))) < sqrt(.Machine$double.eps))

### yeah
m1 <- gamboost(y ~ (bols(x1, intercept = FALSE, df = 1) %+%
                    bols(x2, intercept = FALSE, df = 1)) %X% bols(f, df = 4) +
                    bbs(x1) + bspatial(x1, x2))


### test bmono with categorical covariates
set.seed(781)
x <- as.ordered(gl(4, 50))
y <- rnorm(200, mean = rep(c(1,3,2,4), each = 50), sd = 1)
#plot(x, y)
mod1 <- gamboost(y ~ bols(x))
diff(coef(mod1)[[1]])
mod2 <- gamboost(y ~ bmono(x, lambda2 = 10^15))
stopifnot(abs(coef(mod1)[[1]] - coef(mod2)[[1]])[c(1,4)] < 1e-5)
stopifnot(all(diff(coef(mod2)[[1]]) > - sqrt(.Machine$double.eps)))
