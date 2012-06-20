
require("mboostDevel")

set.seed(290875)

### for boosting hat matrix checks
fm <- Gaussian()
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
data("bodyfat", package = "mboostDevel")
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
stopifnot(max(abs(fitted(mod1) - fitted(mod2))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp))) < sqrt(.Machine$double.eps))

fm1 <- y ~ bbs(x1, df = 3) + bbs(x2, df = 3)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "bbs", dfbase = 3)
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < sqrt(.Machine$double.eps))

fm1 <- y ~ bols(x1) + bols(x2)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "bols")
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < sqrt(.Machine$double.eps))

fm1 <- y ~ btree(x1) + btree(x2)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm2, data = tmp, base = "btree")
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < sqrt(.Machine$double.eps))

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

### predictions:
data("bodyfat", package = "mboostDevel")
amod <- gamboost(DEXfat ~ hipcirc + anthro3a, data = bodyfat, baselearner = "bbs")

agg <- c("none", "sum", "cumsum")
whi <- list(NULL, 1, 2, c(1,2))
for (i in 1:4){
    pred <- vector("list", length=3)
    for (j in 1:3){
        pred[[j]] <- predict(amod, aggregate=agg[j], which = whi[[i]])
    }
    if (i == 1){
        stopifnot(max(abs(pred[[2]] - pred[[3]][,ncol(pred[[3]])]))  < sqrt(.Machine$double.eps))
        if ((pred[[2]] - rowSums(pred[[1]]))[1] - attr(coef(amod), "offset") < sqrt(.Machine$double.eps))
            warning(sQuote("aggregate = sum"), " adds the offset, ", sQuote("aggregate = none"), " doesn't.")
        stopifnot(max(abs(pred[[2]] - rowSums(pred[[1]]) - attr(coef(amod), "offset")))   < sqrt(.Machine$double.eps))
    } else {
        stopifnot(max(abs(pred[[2]] - sapply(pred[[3]], function(obj) obj[,ncol(obj)])))  < sqrt(.Machine$double.eps))
        stopifnot(max(abs(pred[[2]] - sapply(pred[[1]], function(obj) rowSums(obj))))  < sqrt(.Machine$double.eps))
    }
}

# use names in which instead of numbers
agg <- c("none", "sum", "cumsum")
whi <- list(NULL, "hipcirc", "anthro3a", c("hipcirc", "anthro3a"))
for (i in 1:4){
    pred <- vector("list", length=3)
    for (j in 1:3){
        pred[[j]] <- predict(amod, aggregate=agg[j], which = whi[[i]])
    }
    if (i == 1){
        stopifnot(max(abs(pred[[2]] - pred[[3]][,ncol(pred[[3]])]))  < sqrt(.Machine$double.eps))
        if ((pred[[2]] - rowSums(pred[[1]]))[1] - attr(coef(amod), "offset") < sqrt(.Machine$double.eps))
            warning(sQuote("aggregate = sum"), " adds the offset, ", sQuote("aggregate = none"), " doesn't.")
        stopifnot(max(abs(pred[[2]] - rowSums(pred[[1]]) - attr(coef(amod), "offset")))   < sqrt(.Machine$double.eps))
    } else {
        stopifnot(max(abs(pred[[2]] - sapply(pred[[3]], function(obj) obj[,ncol(obj)])))  < sqrt(.Machine$double.eps))
        stopifnot(max(abs(pred[[2]] - sapply(pred[[1]], function(obj) rowSums(obj))))  < sqrt(.Machine$double.eps))
    }
}

# use which to extract all effects for one covariate e.g. for plotting purposes
amod <- gamboost(DEXfat ~ bols(hipcirc, intercept=FALSE) + bbs(hipcirc, df = 1, center=TRUE), data = bodyfat)
stopifnot(ncol(predict(amod, which="hip")) == 2 && all(rowSums(predict(amod, which="hip")) + attr(coef(amod), "offset") - predict(amod) < sqrt(.Machine$double.eps)))

amod <- gamboost(DEXfat ~ hipcirc + anthro3a + kneebreadth,
                 data = bodyfat, baselearner = "bbs")
pr1 <- predict(amod, aggre = "sum", which= 1:2)
foo <- bodyfat[,names(bodyfat) %in% c("hipcirc", "anthro3a", "kneebreadth")]
foo$kneebreadth <- mean(bodyfat$kneebreadth)
pr2 <- predict(amod, aggre = "sum", newdata=foo)
stopifnot(all(diff(rowSums(pr1) - pr2) < sqrt(.Machine$double.eps))) # changes in level are ok
newData <- as.data.frame(rbind(colMeans(bodyfat)[-2], colMeans(bodyfat)[-2]+1*sapply(bodyfat, sd)[-2]))
if (!is.list(pr <- predict(amod, newdata=newData, which=1:2)))
    warning("predict(amod, newdata=newData, which=1:2) does not return a list") # no list but a matrix is returned!
stopifnot(is.list(pr <- predict(amod, newdata=newData, aggregate="cumsum", which=1:2)))
amod[10]
pr <- predict(amod, which=1:3)
stopifnot(ncol(pr) == 3 || all(pr[,ncol] == 0))
amod[100]

# check type argument
set.seed(1907)
x1 <- rnorm(100)
p <- 1/(1 + exp(- 3 * x1))
y <- as.factor(runif(100) < p)
DF <- data.frame(y = y, x1 = x1)

logitBoost <- gamboost(y ~ x1, family = Binomial(),
                 data = DF, baselearner = "bols", control=boost_control(mstop=5000))
logit <- glm(y ~ x1, data=DF, family=binomial)
stopifnot(coef(logitBoost)[[1]][2]*2 - coef(logit)[2]  < sqrt(.Machine$double.eps)) # * 2 as we use y = {-1, 1}

pr <- predict(logitBoost)
pr2 <- predict(logit)
stopifnot(pr * 2 - pr2 < 1e-5)  # * 2 as we use y = {-1, 1}

pr <- predict(logitBoost, type="class")
pr2 <- predict(logit, type="response") > 0.5
foo <- table(pr, pr2)
stopifnot(foo[1,2] + foo[2,1] == 0)

pr <- predict(logitBoost, type="response")
pr2 <- predict(logit, type="response")
stopifnot(pr - pr2  < sqrt(.Machine$double.eps))

### coefficients:
data("bodyfat", package = "mboostDevel")
amod <- gamboost(DEXfat ~ hipcirc + anthro3a + kneebreadth,
                 data = bodyfat, baselearner = "bbs")
stopifnot(length(coef(amod)) == 3)
amod[10]
stopifnot(length(coef(amod)) == 2)
stopifnot(length(coef(amod, which=1:3)) == 3)

### cyclic covariates
x <- seq(from = 0, to = 2*pi, length = 100)
y <- sin(x) + rnorm(length(x), sd = 0.5)
mod <- gamboost(y ~ bbs(x, cyclic = TRUE))
stopifnot(diff(fitted(mod)[c(1, 100)]) == 0)

### buser
set.seed(1907)
x <- rnorm(100)
y <- rnorm(100, mean = x^2, sd = 0.1)
mod1 <- gamboost(y ~ bbs(x))
X <- extract(bbs(x))
K <- extract(bbs(x), "penalty")
mod2 <- gamboost(y ~ buser(X, K))
stopifnot(max(abs(predict(mod1) - predict(mod2))) < sqrt(.Machine$double.eps))

z <- sample(1:2, 100, replace=TRUE)
y[z == 2] <- rnorm(100, mean = - x^2, sd = 0.1)[z == 2]
z <- as.factor(z)
mod3 <- gamboost(y ~ bbs(x) + bbs(x, by = z),
                 control = boost_control(mstop = 1000))
X <- extract(bbs(x))
K <- extract(bbs(x), "penalty")
mod4 <- gamboost(y ~  buser(X, K) + buser(X, K, by = z),
                 control = boost_control(mstop = 1000))
stopifnot(max(abs(predict(mod3) - predict(mod4))) < sqrt(.Machine$double.eps))

y <- rnorm(100, mean = as.numeric(z), sd = 0.1)
mod5 <- gamboost(y ~ bols(z))
X <- extract(bols(z))
K <- extract(bols(z), "penalty")
index <- extract(bols(z), "index")
mod6 <- gamboost(y ~  buser(X, K, lambda = 0, index = index))
mod6a <- gamboost(y ~  buser(X, index = index))
stopifnot(max(abs(predict(mod5) - predict(mod6))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(predict(mod5) - predict(mod6a))) < sqrt(.Machine$double.eps))

z <- sample(1:3, 100, replace = TRUE)
y <- rnorm(100, mean = z, sd = 0.1)
z <- as.ordered(z)
mod7 <- gamboost(y ~ bols(z))
X <- extract(bols(z))
K <- extract(bols(z), "penalty")
index <- extract(bols(z), "index")
mod8 <- gamboost(y ~  buser(X, K, lambda = 0, index = index))
stopifnot(max(abs(predict(mod7) - predict(mod8))) < sqrt(.Machine$double.eps))

y[z == 1] <- rnorm(100, mean = x^2, sd = 0.1)[z == 1]
y[z == 2] <- rnorm(100, mean = - x^2, sd = 0.1)[z == 2]
y[z == 3] <- rnorm(100, mean = x, sd = 0.1)[z == 3]
mod9 <- gamboost(y ~ bbs(x) + bbs(x, by = z),
                 control = boost_control(mstop = 1000))
X <- extract(bbs(x))
K <- extract(bbs(x), "penalty")
mod10 <- gamboost(y ~  buser(X, K) + buser(X, K, by = z),
                 control = boost_control(mstop = 1000))
stopifnot(max(abs(predict(mod9) - predict(mod10))) < sqrt(.Machine$double.eps))
