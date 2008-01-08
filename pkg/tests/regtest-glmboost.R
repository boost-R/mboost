
require("mboost")

set.seed(290875)

dgp <- function(n = 100, beta = rep(0, 10), sd = 1) {
    p <- length(beta) - 1
    x <- cbind(1, matrix(runif(n * p), ncol = p))
    lp <- x %*% beta
    y <- lp + rnorm(n, sd = sd)
    ls <- data.frame(y = y, x[,-1])
    attr(ls, "lp") <- lp
    ls
}

### well-defined problem
mydf <- dgp(beta = c(1, 2.5, rep(0, 2)))

### for easy comparison with lm
fm <- GaussReg()
fm@offset <- function(y, w) 0

mydf.gb <- glmboost(y ~ ., data = mydf, family = fm, 
                    control = boost_control(mstop = 1000, nu = 1))
mydf.lm <- lm(y ~ ., data = mydf)

### compare coefficients
stopifnot(max(abs(coef(mydf.gb) - coef(mydf.lm))) < 1e-10)

### a little bit more difficult
mydf <- dgp(beta = c(1, 2.5, rep(0, 38)))

mydf.gb <- glmboost(y ~ ., data = mydf, family = fm, 
                    control = boost_control(mstop = 1000, nu = 1))
aic <- AIC(mydf.gb, method = "corrected")
ht <- hatvalues(mydf.gb)
mstop(aic)
mydf.lm <- lm(y ~ ., data = mydf)

### compare coefficients
which(abs(coef(mydf.lm)) < abs(coef(mydf.gb[mstop(aic)])))

#### check boosting hat matrix and subsetting / predict
stopifnot(isTRUE(all.equal(drop(attr(ht, "hatmatrix") %*% mydf$y),
                           as.vector(predict(mydf.gb)))))
ht25 <- hatvalues(mydf.gb[25])
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% mydf$y),
                           as.vector(predict(mydf.gb[25])))))
stopifnot(isTRUE(all.equal(drop(attr(ht25, "hatmatrix") %*% mydf$y),
                           as.vector(fitted(mydf.gb[25])))))

### a simple two-dimensional example from `glmboost.Rd'
data("cars")
cars.gb <- glmboost(dist ~ speed, data = cars, family = fm, 
                    control = boost_control(mstop = 1000, nu = 1))
cars.gb

### coefficients should coincide
cf <- coef(cars.gb)
attr(cf, "offset") <- NULL
stopifnot(all.equal(cf, coef(lm(dist ~ speed, data = cars))))

### logistic regression
mydf <- data.frame(x = runif(100), z = rnorm(100), 
                   y = factor(c(rep(0, 30), rep(1, 70))))
bmod <- glmboost(y ~ x + z, data = mydf, family = Binomial(), 
                 control = boost_control(mstop = 1000, nu = 1))
gmod <- glm(y ~ x + z, data = mydf, family = binomial())
llg <- logLik(gmod)
attributes(llg) <- NULL
stopifnot(all.equal(logLik(bmod), llg))
stopifnot(max(abs(predict(gmod, type = "link")/2 - fitted(bmod))) < 
                  sqrt(.Machine$double.eps))
cfb <- (coef(bmod) + c(bmod$offset, 0, 0)) * 2
attr(cfb, "offset") <- NULL
stopifnot(all.equal(cfb, coef(gmod)))
aic <- AIC(bmod, "classical")
stopifnot(abs(AIC(gmod) - attr(aic, "AIC")[mstop(bmod)]) < 1e-5)

### weighted least squares problem

x <- runif(100)  
df <- data.frame(y = 2 + 3 * x + rnorm(length(x)),
                 x = x, z = runif(length(x)),
                 w = runif(length(x)) * 10)

### linear model, classical fit
lmmod <- lm(y ~ x + z, data = df, weights = w)

### linear model, boosting fit
lmb <- glmboost(y ~ x + z, data = df, weights = df$w,
                control = boost_control(mstop = 5000, nu = 1))

### compare fitted values
stopifnot(max(abs(fitted(lmmod) -fitted(lmb))) < sqrt(.Machine$double.eps))

### compare hat matrices
stopifnot(max(abs(hatvalues(lmmod) - hatvalues(lmb))) < sqrt(.Machine$double.eps))      

### compare boosting hat matrix with fitted values
stopifnot(max(abs(attr(hatvalues(lmb), "hatmatrix") %*% (df$y - lmb$offset) + lmb$offset -
        fitted(lmb))) < sqrt(.Machine$double.eps)) 

### Cox model (check for CoxPH family)
if (require("survival")) {

    test <- data.frame(time = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                       event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                       x     = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0))

    stopifnot(all.equal(coef(cx <- coxph(Surv(time, event) ~ x, data = test, method = "breslow")),
                       coef(gl <- glmboost(Surv(time, event) ~ x, data = test,
                       family = CoxPH(), 
                       control = boost_control(mstop = 2000, nu = 1)))[2]))

    stopifnot(all.equal(cx$loglik[2], logLik(gl)))

    indx <- c(1, 1, 1, 2:10)
    w <- tabulate(indx)

    stopifnot(all.equal(coef(cx <- coxph(Surv(time, event) ~ x, data = test, weights = w, 
                                   method = "breslow")),
                       coef(gl <- glmboost(Surv(time, event) ~ x, data = test, weights = w,
                       family = CoxPH(), 
                       control = boost_control(mstop = 200, nu = 1)))[2]))

    stopifnot(all.equal(cx$loglik[2], logLik(gl)))

    indx <- c(1, 1, 1, 3:10)
    w <- tabulate(indx)

    stopifnot(all.equal(coef(cx <- coxph(Surv(time, event) ~ x, data = test[indx,], 
                                   method = "breslow")),
                       coef(gl <- glmboost(Surv(time, event) ~ x, data = test, weights = w,
                       family = CoxPH(), 
                       control = boost_control(mstop = 1000)))[2]))

    stopifnot(all.equal(cx$loglik[2], logLik(gl)))
}

### check centering
y <- rnorm(20)
xn <- rnorm(20)
xnm <- xn - mean(xn)
xf <- gl(2, 10)
gc <- glmboost(y ~ xn + xf, control = boost_control(center = TRUE))
g <- glmboost(y ~ xnm + xf)
cgc <- coef(gc)
cg <- coef(g)
names(cgc) <- NULL
names(cg) <- NULL
stopifnot(all.equal(cgc, cg))
stopifnot(all.equal(mstop(AIC(gc, "corrected")), mstop(AIC(g, "corrected"))))

pc1 <- predict(gc)
pc2 <- predict(gc, newdata = data.frame(xn = xn, xf = xf))
pc3 <- predict(g)
stopifnot(all.equal(pc1, pc2))
stopifnot(all.equal(pc2, pc3))

### Poisson models

df <- data.frame(x1 = runif(100), x2 = runif(100))
f <- -1 + 3 * df$x2
df$y <- round(exp(f) )
ctrl <- boost_control(mstop = 2000, nu = 0.1)

gmod <- glm(y ~ x1 + x2, data = df, family = poisson())
gbmod <- glmboost(y ~ x1 + x2, data = df, family = Poisson(), control = ctrl)

llg <- logLik(gmod)
attributes(llg) <- NULL
stopifnot(all.equal(logLik(gbmod), llg))

### hat matrix is only approximate!
stopifnot(abs(AIC(gmod) - attr(AIC(gbmod, "classical"), "AIC")[mstop(gbmod)]) < 1)

stopifnot(max(abs(predict(gmod) -  predict(gbmod))) < 1e-4)
