
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
                           as.vector(predict(cars.gb)))))
ht25 <- hatvalues(cars.gb[25])
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% cars$dist),
                           as.vector(predict(cars.gb[25])))))
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% cars$dist),
                           as.vector(fitted(cars.gb[25])))))

### check boosting hat matrix with multiple independent variables
### and weights
data("bodyfat", package = "mboost")
bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4 - 1
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
ga <- gamboost(y ~ x.2 - 1, data = df,
               control = boost_control(mstop = 100, nu = 1))
stopin(fitted(lm(y ~ x.2 - 1, data = df)), fitted(ga))

### univariate model involving sin transformation
df <- data.frame(y = sin(x[,1]), x = x)
ga <- gamboost(y ~ x.1 - 1, data = df, 
               control = boost_control(mstop = 100, nu = 1))
stopin(fitted(lm(y ~ sin(x.1) - 1, data = df)), fitted(ga))

### bivariate model: linear and sin
df <- data.frame(y = sin(x[,1]) + 3*x[,2], x = x)
ga <- gamboost(y ~ x.1 + x.2 - 1, data = df, 
               control = boost_control(mstop = 100, nu = 1))
stopin(fitted(lm(y ~ sin(x.1) + x.2 - 1, data = df)), fitted(ga))
ga <- gamboost(y ~ x.1 + x.2 - 1, data = df, dfbase = c(4, 1), 
               control = boost_control(mstop = 100, nu = 1))
stopin(fitted(lm(y ~ sin(x.1) + x.2 - 1, data = df)), fitted(ga))

### ANCOVA model
df <- data.frame(y = 3 * x[,2] + (1:4)[xf], x = x)
ga <- gamboost(y ~ xf + x.2 - 1, data = df, 
               control = boost_control(mstop = 100, nu = 1))
stopin(fitted(lm(y ~ xf + x.2 - 1, data = df)), fitted(ga))
ga <- gamboost(y ~ xf + sin(x.1) + x.2, data = df, 
               dfbase = c(1, 1, 4, 1),
               control = boost_control(mstop = 100, nu = 1))
stopin(fitted(lm(y ~ xf + sin(x.1) + x.2, data = df)), fitted(ga))


### check centering
y <- rnorm(20)
xn <- rnorm(20)
xnm <- xn - mean(xn)
xf <- gl(2, 10)
gc <- gamboost(y ~ xn + xf, control = boost_control(center = TRUE))
g <- gamboost(y ~ xnm + xf)
cgc <- coef(gc)
cg <- coef(g)  
names(cgc) <- NULL
names(cg) <- NULL 
stopifnot(all.equal(cgc, cg))

pc1 <- predict(gc)
pc2 <- predict(gc, newdata = data.frame(xn = xn, xf = xf))
pc3 <- predict(g)
stopifnot(all.equal(pc1, pc2))
stopifnot(all.equal(pc2, pc3))

### formula interfaces
tmp <- data.frame(x1 = runif(100), x2 = runif(100), y = rnorm(100))
fm1 <- y ~ bss(x1, df = 3) + bss(x2, df = 3)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm1, data = tmp, base = "bss", dfbase = 3)
stopifnot(max(abs(fitted(mod1) - fitted(mod2))) < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp))) < .Machine$double.eps)

fm1 <- y ~ bbs(x1, df = 3) + bbs(x2, df = 3)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm1, data = tmp, base = "bbs", dfbase = 3)
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < .Machine$double.eps)

fm1 <- y ~ bols(x1) + bols(x2)
fm2 <- y ~ x1 + x2
mod1 <- gamboost(fm1, data = tmp)
mod2 <- gamboost(fm1, data = tmp, base = "bols")
stopifnot(max(abs(fitted(mod1) - fitted(mod2)))  < .Machine$double.eps)
stopifnot(max(abs(predict(mod1, newdata = tmp) - predict(mod2, newdata = tmp)))  < .Machine$double.eps)


### Weibull model with scale parameter estimation

### random numbers from extreme value distribution
rextrval <- function(x) log( -log(1-x) )

n <- 300
sigma <- 0.5
u <- runif(n)
u.1 <- runif(n)
w <- rextrval(u)
w.1 <- rextrval(u.1)

x1 <- runif(n,-3,3)
x1.1 <- runif(n,-3,3)
x2 <- x1 + runif(n,-3,3)
x2.1 <- x1.1 + runif(n,-3,3)
survtime <- exp(sin(x1) + 0.5*cos(x2) + sigma*w)
censtime <- exp(sin(x1.1) + 0.5*cos(x2.1) + sigma*w.1)
event <- survtime<censtime
stime <- pmin(survtime,censtime)

###
ctrl <- boost_control(center=T,mstop=1000)
model1 <- gamboost(Surv(stime,event)~bbs(x1, df=3, knots=40)+bbs(x2, df=3,
    knots=40), family=Weib(), control=ctrl)
par(mfrow=c(1,2))
plot(model1, ask=F)
model1 <- gamboost(Surv(stime,event)~bss(x1, df=3)+bss(x2, df=3),
family=Weib(), control=ctrl)
par(mfrow=c(1,2))
plot(model1, ask=F)


### Log logistic model with scale parameter estimation

sigma <- 0.5
n <- 300
w <- rlogis(n)
w.1 <- rlogis(n)

x1 <- runif(n,-3,3)
x1.1 <- runif(n,-3,3)
x2 <- x1 + runif(n,-3,3)
x2.1 <- x1.1 + runif(n,-3,3)
survtime <- exp(sin(x1) + 0.5*cos(x2) + sigma*w)
censtime <- exp(sin(x1.1) + 0.5*cos(x2.1) + sigma*w.1)
event <- survtime<censtime
stime <- pmin(survtime,censtime)

###
ctrl <- boost_control(center=T,mstop=1000)
model1 <- gamboost(Surv(stime,event)~bbs(x1, df=3, knots=40)+bbs(x2, df=3,
    knots=40), family=Loglog(), control=ctrl)
par(mfrow=c(1,2))
plot(model1, ask=F)
model1 <- gamboost(Surv(stime,event)~bss(x1, df=3)+bss(x2, df=3),
family=Loglog(), control=ctrl)
par(mfrow=c(1,2))
plot(model1, ask=F)

### Log normal model with scale parameter estimation

sigma <- 0.5
n <- 300
w <- rnorm(n)
w.1 <- rnorm(n)

x1 <- runif(n,-3,3)
x1.1 <- runif(n,-3,3)
x2 <- x1 + runif(n,-3,3)
x2.1 <- x1.1 + runif(n,-3,3)
survtime <- exp(sin(x1) + 0.5*cos(x2) + sigma*w)
censtime <- exp(sin(x1.1) + 0.5*cos(x2.1) + sigma*w.1)
event <- survtime<censtime
stime <- pmin(survtime,censtime)

###
ctrl <- boost_control(center=T,mstop=1000)
model1 <- gamboost(Surv(stime,event)~bbs(x1, df=3, knots=40)+bbs(x2, df=3,
    knots=40), family=LogNormal(), control=ctrl)
par(mfrow=c(1,2))
plot(model1, ask=F)
model1 <- gamboost(Surv(stime,event)~bss(x1, df=3)+bss(x2, df=3),
family=LogNormal(), control=ctrl)
par(mfrow=c(1,2))
plot(model1, ask=F)

