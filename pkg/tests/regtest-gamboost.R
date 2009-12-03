
require("mboost")

set.seed(290875)

### for boosting hat matrix checks
fm <- GaussReg()
fm@offset <- function(y, w) 0

### a simple two-dimensional example from `gamboost.Rd'
data("cars")
cars.gb <- gamboost(dist ~ speed, data = cars, df = 4, family = fm,
                    control = boost_control(mstop = 50))
cars.gb
aic <- AIC(cars.gb, method = "corrected")
aic

ht <- hatvalues(cars.gb)

### plot fit
plot(dist ~ speed, data = cars)
lines(cars$speed, predict(cars.gb[mstop(AIC(cars.gb))]), col = "red")
lines(cars$speed, predict(smooth.spline(cars$speed, cars$dist), cars$speed)$y, 
      col = "green")

#### check boosting hat matrix and subsetting / predict
stopifnot(isTRUE(all.equal(drop(attr(ht, "hatmatrix") %*% cars$dist),
                           as.vector(predict(cars.gb[50])))))
ht25 <- hatvalues(cars.gb[25])
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% cars$dist),
                           as.vector(predict(cars.gb[25])))))
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% cars$dist),
                           as.vector(fitted(cars.gb[25])))))

### check boosting hat matrix with multiple independent variables
### and weights
data("bodyfat", package = "mboost")
bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4 
indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
bodyfat[indep] <- lapply(bodyfat[indep], function(x) x - mean(x))
bf_gam <- gamboost(bffm, data = bodyfat, control = boost_control(mstop = 10), 
                   weights = runif(nrow(bodyfat)) * 10)
### aic <- AIC(bf_gam)
ht <- hatvalues(bf_gam)

off <- bf_gam$offset
u <- bf_gam$ustart

stopifnot(isTRUE(all.equal(drop(attr(ht, "hatmatrix") %*% u + off),
                           as.vector(predict(bf_gam)))))
stopifnot(isTRUE(all.equal(drop(attr(ht, "hatmatrix") %*% u + off),
                           as.vector(fitted(bf_gam)))))


### compare `gamboost' with `lm' in cases where this is actually possible
set.seed(290875)
x <- matrix(runif(1000) * 10, ncol = 10)
xf <- gl(4, nrow(x)/4)

### OK, we need to allow for some small differences (larger mstop values
### would fix this)
stopin <- function(x, y) stopifnot(max(abs(x - y)) < 0.1)

### univariate linear model
df <- data.frame(y = 3*x[,2], x = x)
ga <- gamboost(y ~ x.2, data = df,
               control = boost_control(mstop = 500, nu = 1))
stopin(fitted(lm(y ~ x.2 - 1, data = df)), fitted(ga))

### univariate model involving sin transformation
df <- data.frame(y = sin(x[,1]), x = x)
ga <- gamboost(y ~ x.1, data = df, 
               control = boost_control(mstop = 500, nu = 1))
stopin(fitted(lm(y ~ sin(x.1) - 1, data = df)), fitted(ga))

### bivariate model: linear and sin
df <- data.frame(y = sin(x[,1]) + 3*x[,2], x = x)
ga <- gamboost(y ~ x.1 + x.2, data = df, 
               control = boost_control(mstop = 500, nu = 1))
stopin(fitted(lm(y ~ sin(x.1) + x.2 - 1, data = df)), fitted(ga))
ga <- gamboost(y ~ x.1 + bols(x.2), data = df, 
               control = boost_control(mstop = 500, nu = 1))
stopin(fitted(lm(y ~ sin(x.1) + x.2 - 1, data = df)), fitted(ga))

### ANCOVA model
df <- data.frame(y = 3 * x[,2] + (1:4)[xf], x = x)
ga <- gamboost(y ~ xf + x.2, data = df, 
               control = boost_control(mstop = 500, nu = 1))
stopin(fitted(lm(y ~ xf + x.2 - 1, data = df)), fitted(ga))
ga <- gamboost(y ~ xf + sin(x.1) + x.2, data = df, 
               dfbase = 1,
               control = boost_control(mstop = 500, nu = 1))
stopin(fitted(lm(y ~ xf + sin(x.1) + x.2, data = df)), fitted(ga))

### check centering
y <- rnorm(200)
xn <- rnorm(200)
xnm <- xn - mean(xn)
xf <- gl(2, 100)
gc <- gamboost(y ~ xn + xf)
g <- gamboost(y ~ xnm + xf)
stopifnot(max(abs(fitted(gc) - fitted(g))) < 1 / 10000)

pc1 <- predict(gc)
pc2 <- predict(gc, newdata = data.frame(xn = xn, xf = xf))
pc3 <- predict(g)
stopifnot(all.equal(pc1, pc2))
stopifnot(max(abs(pc2 - pc3)) < 1 / 10000)

### formula interfaces
tmp <- data.frame(x1 = runif(100), x2 = runif(100), y = rnorm(100))
fm1 <- y ~ bbs(x1, df = 3) + bbs(x2, df = 3)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "bss", dfbase = 3)
stopifnot(max(abs(fitted(mod1) - fitted(mod2))) < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp))) < .Machine$double.eps)

fm1 <- y ~ bbs(x1, df = 3) + bbs(x2, df = 3)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "bbs", dfbase = 3)
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < .Machine$double.eps)

fm1 <- y ~ bols(x1) + bols(x2)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "bols")
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < .Machine$double.eps)

fm1 <- y ~ btree(x1) + btree(x2)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "btree")
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < .Machine$double.eps)

## Cox model

fit2 <- gamboost(Surv(futime, fustat) ~ bbs(age) +
    bols(resid.ds) + bols(rx) + bols(ecog.ps), data = ovarian, 
    family = CoxPH(), control = boost_control(mstop = 1000))

A2 <- survFit(fit2)
A2

newdata <- ovarian[c(1,3,12),]
A2 <- survFit(fit2, newdata = newdata)
A2

### gamboost with explicit intercept
df <- data.frame(x = 1:100, y = rnorm(1:100), int = rep(1, 100))
mod <- gamboost(y ~ bols(int, intercept = FALSE) + bols(x, intercept = FALSE), data = df, 
                control = boost_control(mstop = 2500))
cf <- unlist(coef(mod))
cf[1] <- cf[1] + mod$offset
tmp <- max(abs(cf - coef(lm(y ~ x, data = df))))
stopifnot(tmp < 1e-5)
tmp <- max(abs(fitted(mod) - fitted(lm(y ~ x, data = df))))
stopifnot(tmp < 1e-5)

### predictions <FIXME>: more tests </FIXME>
data("bodyfat", package = "mboost")
amod <- gamboost(DEXfat ~ hipcirc + anthro3a + kneebreadth, 
                 data = bodyfat, baselearner = "bbs")
pr <- predict(amod, aggre = "cumsum", which = 1:2)
pr <- predict(amod, aggre = "none", which = 1:2)
